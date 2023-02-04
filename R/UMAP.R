# work.path <- "i:/genomicdata/Research/BatchEffect";setwd(work.path)
# article.path <- file.path(work.path, "Literatures")
# data.path <- file.path(work.path, "InputData")
# pkg.path <- file.path(work.path, "Packages")
# code.path <- file.path(work.path, "Codes")
# res.path <- file.path(work.path, "Results")
# fig.path <- file.path(work.path, "Figures")
#
# res.tmp.path <- file.path(res.path, "Visualization")
# fig.tmp.path <- file.path(fig.path, "Visualization")
#
# invisible(lapply(ls()[grep("path", ls())], function(x){
#   if (!dir.exists(get(x))) dir.create(get(x))
# }))

# BiocManager::install("densvis")

# library(umap)
# library(densvis)
# library(Seurat)
# library(SeuratObject)
# library(leiden)

# edit(umap)
# edit(RunUMAP)
# methods(RunUMAP)
# View(Seurat:::RunUMAP.default)
# View(Seurat:::RunUMAP.Graph)
# View(Seurat:::FindNeighbors.default)
# View(Seurat:::ComputeSNN)
# View(Seurat:::FindClusters.default)
# View(Seurat:::RunModularityClustering)
# View(leiden:::leiden.matrix)

# A test ------------------------------------------------------------------


# seu <- readRDS("G:/singlecell/pbmc.rds")
# seu <- FindNeighbors(seu)
# seu <- FindClusters(seu, resolution = 0.4)
# seu <- RunUMAP(object = seu, graph = "RNA_snn", umap.method = "umap-learn")
# seu <- RunUMAP(object = seu, reduction = "pca", dims = )
# neigh <- seu@graphs
# neighborlist <- list("idx" = Indices(neigh),
#                      "dist" = Distances(neigh))


# UMAP R code -------------------------------------------------------------

# iris.umap = umap(d = iris[,1:4],
#                  config = umap.defaults,
#                  method = "naive")
# plot(iris.umap$layout)
#
# d = pcs
# method = "naive"
# config = umap.defaults
#
# method = config$method = match.arg(arg = method, choices = c("naive", "umap-learn"))
# config = umap:::umap.check.config(config)
# d = umap:::umap.prep.input(d, config)
# old.seed = umap:::get.global.seed()
#
# # implementations = c(naive = umap.naive, `umap-learn` = umap.learn)
# # result = implementations[[method]](d, config)
# # View(umap:::umap.naive)
#
# if (is.na(config$a) | is.na(config$b))
#   config[c("a", "b")] = umap:::find.ab.params(config$spread, config$min_dist)
# if (is.na(config$random_state))
#   config$random_state = as.integer(runif(1, 0, 2^30))
# if (is(d, "Matrix"))
#   d = Matrix::as.matrix(d)
#
# knn = NULL
# if ("knn" %in% names(config) & is(config$knn, "umap.knn"))
#   knn = config$knn
# if (is.null(knn))
#   knn = umap:::knn.info(d, config)
# graph = umap:::naive.fuzzy.simplicial.set(knn, config)
# embedding = umap:::make.initial.embedding(graph$n.elements, config, graph)
# embedding = umap:::naive.simplicial.set.embedding(graph, embedding, config)
# embedding = umap:::center.embedding(embedding)
# list(layout = embedding, data = d, knn = knn, config = config)


# Try to make UMAP based on Shared Nearest Neighbor ----------------------


# Rcpp::sourceCpp("/share/genomics/cwx/Research/snnMAP/Codes/snnMAPhelper.cpp")
# Rcpp::sourceCpp("/share/genomics/cwx/Research/snnMAP/Codes/iterator.cpp")

# library(parallel)

RunSnnMAP <- function(seu, dims, mat = NULL, k.param = 20, complexity = 2,
                      ncores = getOption("mc.cores", 2L), n.iter = 10){
  # mat: n cell * p feature
  if (is.null(mat)){
    pcs <- seu@reductions$pca@cell.embeddings[, dims]
  }else{
    pcs <- mat
  }
  snn <- SeuratObject::as.sparse(Seurat:::FindNeighbors(pcs, k.param = k.param)$snn)
  snn <- as(snn, "dgTMatrix")

  method_args = list("gamma" = 1, "approx_pow" = F)
  # method_args[c("a", "b")] = uwot:::find_ab_params(config$spread, config$min_dist)
  n_epochs = 500; n_components = 2; n_vertices = nrow(pcs); negative_sample_rate = 5; grain_size = 1; alpha = 1
  full_opt_args <- uwot:::get_opt_args(opt_args = NULL, alpha)
  epoch_callback = NULL; pcg_rand = T; batch = F; n_sgd_threads = ncores

  message("calculate the graph for ", nrow(pcs), " cells ")
  snn.df <- data.frame("from" = snn@i+1, "to" = snn@j+1, "snn" = snn@x)
  snn.df <- dplyr::arrange(snn.df, snn.df$from, snn.df$to)
  snn.df$dist <- rcpp_distance(mat = pcs, X = snn.df$from, Y = snn.df$to)
  snn.df <- split(snn.df, snn.df$from)

  ## set alternative model
  formula.list <- lapply(0:complexity, function(com){
    if (com==0){
      return("snn ~ exp(-b0*dist)")
    }else{
      paste0("snn ~ exp(-b0*dist) + ", paste(c(sapply(1:com, function(i){
        paste0("a", i, "*exp(b", i, "*dist)")
      })), collapse = " + "))
    }
  })

  options("mc.cores" = ncores)
  message("Estimating local density for each cell\n")
  fit <- mclapply(snn.df, function(data){
    # cat(unique(data$from), "\t")
    data$score <- data$snn * data$dist
    data <- dplyr::arrange(data, -data$score)

    ## Take the rough parameter estimation as start
    param.list <- lapply(0:complexity, function(com){
      param <- list()
      param[["b0"]] <- log(1/mean(data$snn))/mean(data$dist)
      if (com>0){
        for (i in 1:com){
          param[[paste0("a", i)]] = 0.01
          param[[paste0("b", i)]] = log(data$snn[i]/0.01)/data$dist[i]
        }
      }
      return(param)
    })

    ## Estimate the parameter for each model
    model.list <- mapply(function(f, st){
      tryCatch({gslnls::gsl_nls(fn = formula(f), data = data, start = unlist(st), control = list(maxiter = n.iter))},
      # tryCatch({nls(formula = f, data = data, start = st)},
               error = function(e){})
    }, formula.list, param.list, SIMPLIFY = F)
    model.list <- model.list[!unlist(lapply(model.list, is.null))]

    ## Avaerage the estimation
    b0.coef <- unlist(lapply(model.list, function(x){coefficients(x)["b0"]}))
    b0.weight <- unlist(lapply(model.list, function(x){-0.5*AIC(x)}))
    b0.weight[b0.weight<0] = 0
    b0.weight <- b0.weight-min(b0.weight)
    b0.weight <- exp(b0.weight)
    b0.weight[b0.coef<0]=0
    return(c("b0" = sum(b0.coef*b0.weight/sum(b0.weight))))
  })
  message("Fail to estimate the density ", sum(unlist(lapply(fit, is.na))), " cell")

  ## Estimate parameters of the low dimensional distribution
  param <- mcmapply(function(x, y){
    x0 = median(x$dist); leb0 = log(exp(y["b0"]*x0)-1)
    a0 = min(exp(y["b0"])-1, exp(y["b0"]*x0)-1)
    re = list("a" = a0,  "b" = (leb0-log(a0))/2/log(x0))
    fe = tryCatch({coefficients(
      nls(formula = yv ~ 1/(1 + a * xv^(2 * b)) - exp(-y["b0"]*xv),
          data = data.frame("xv" = x$dist, "yv" = 0),
          start = re)
    )}, error = function(e){})
    setNames(object = if (is.null(fe)) unlist(re) else fe,
             nm = c("a", "b"))
  }, snn.df, fit, SIMPLIFY = F)

  message("Fail to estimate the a and b for ", length(which(unlist(lapply(param, length))==0)), " cells")
  method_args$ai <- unlist(lapply(param, function(x) x["a"]))
  method_args$bi <- unlist(lapply(param, function(x) x["b"]))
  method_args$ndim <- n_components

  message("calculate embeddings\n")

  ## get snnMAP embedding
  embedding = uwot:::spectral_init(snn, ndim = n_components)
  embedding <- t(embedding)
  input <- list(head_embedding = embedding, tail_embedding = NULL,
                positive_head = snn@i, positive_tail = snn@j, positive_ptr = NA,
                n_epochs = n_epochs, epochs_per_sample = uwot:::make_epochs_per_sample(snn@x, n_epochs),
                n_head_vertices = n_vertices, n_tail_vertices = n_vertices,
                method = "snnmap", method_args = method_args, initial_alpha = alpha,
                opt_args = full_opt_args, negative_sample_rate = negative_sample_rate,
                pcg_rand = pcg_rand, batch = batch, n_threads = n_sgd_threads,
                grain_size = grain_size, move_other = TRUE, epoch_callback = epoch_callback,
                verbose = T)
  # str(input$method_args);summary(input$method_args$ai);summary(input$method_args$bi)
  embedding = do.call(optimize_layout_r, input)
  embedding <- t(embedding)
  # plot(embedding, cex = .1, pch = 19)
  rownames(embedding) = rownames(pcs)
  colnames(embedding) = c("snnMAP_1", "snnMAP_2")

  if (is.null(mat)){
    seu@reductions$snnMAP <- Seurat::CreateDimReducObject(
      embeddings = embedding,
      assay = "RNA",
      global = T,
      key = "snnMAP_"
    )
    return(seu)
  }else{
    return(embedding)
  }
}


# Python Code --------------------------------------------------------------------


# reticulate::use_python("D:/Software/python39/python.exe")
# reticulate::py_config()
# Sys.setenv("NUMBA_DISABLE_INTEL_SVML" = 1)
#
# library(reticulate)
# np <- import("numpy", delay_load = TRUE)
# sp <- import("scipy", delay_load = TRUE)
# sklearn <- import("sklearn", delay_load = TRUE)
# umap <- import("umap", delay_load = TRUE)
#
# object <- seu
# diag(x = object) <- 0
# data <- object
# object <- sp$sparse$coo_matrix(arg1 = object)
# ab.params <- umap$umap_$find_ab_params(spread = spread,
#                                        min_dist = min.dist)
# a <- a %||% ab.params[[1]]
# b <- b %||% ab.params[[2]]
# n.epochs <- n.epochs %||% 0L
# random.state <- sklearn$utils$check_random_state(seed = as.integer(x = seed.use))
# umap.args <- list(data = data, graph = object, n_components = n.components,
#                   initial_alpha = learning.rate, a = a, b = b, gamma = repulsion.strength,
#                   negative_sample_rate = negative.sample.rate, n_epochs = as.integer(x = n.epochs),
#                   random_state = random.state, init = "spectral", metric = metric,
#                   metric_kwds = metric.kwds, verbose = verbose)
# if (numeric_version(x = umap$pkg_resources$get_distribution("umap-learn")$version) >=
#     numeric_version(x = "0.5.0")) {
#   umap.args <- c(umap.args, list(densmap = densmap, densmap_kwds = densmap.kwds,
#                                  output_dens = FALSE))
# }
# embeddings <- do.call(what = umap$umap_$simplicial_set_embedding,
#                       args = umap.args)
# if (length(x = embeddings) == 2) {
#   embeddings <- embeddings[[1]]
# }
# rownames(x = embeddings) <- colnames(x = data)
# colnames(x = embeddings) <- paste0("UMAP_", 1 : n.components)
# embeddings <- scale(x = embeddings, scale = FALSE)
# umap <- CreateDimReducObject(embeddings = embeddings, key = reduction.key,
#                              assay = assay, global = TRUE)


# denSMAP -----------------------------------------------------------------


RunDenSMAP <- function(seu, mat = NULL, dims, verbose = F){
  umap.import <- reticulate::import("umap")
  umap.args <- list(n_neighbors = 30L, n_components = 2L, n_epochs = 750L,
                    metric = "euclidean", learning_rate = 1,
                    min_dist = 0.1, spread = 1, set_op_mix_ratio = 1,
                    local_connectivity = 1L, repulsion_strength = 1,
                    negative_sample_rate = 5L, transform_queue_size = 4, random_state = NULL,
                    a = NULL, b = NULL, metric_kwds = NULL, angular_rp_forest = F, target_n_neighbors = -1,
                    target_weight = 0.5, disconnection_distance = NULL,
                    verbose = verbose,
                    densmap = T,
                    dens_lambda = 0.1, dens_frac = 0.3,
                    dens_var_shift = 0.1, output_dens = T)
  umap <- do.call(what = umap.import$UMAP, args = umap.args)

  if (is.null(mat)){
    pcs = as.matrix(x = seu@reductions$pca@cell.embeddings[, paste0("PC_", dims)])
  }else{
    pcs = mat
  }


  embedding <- umap$fit_transform(pcs)
  # embedding <- densvis::densmap(seu@reductions$pca@cell.embeddings[, paste0("PC_", dims)])
  colnames(embedding[[1]]) <- c("denSMAP_1", "densMAP_2")
  rownames(embedding[[1]]) <- rownames(pcs)

  if (is.null(mat)){
    seu@reductions$denSMAP <- Seurat::CreateDimReducObject(
      embeddings = embedding[[1]],
      assay = "RNA",
      global = T,
      key = "denSMAP_"
    )
    seu@reductions$denSMAP@misc <- embedding[c(2,3)]
    return(seu)
  }else{
    return(embedding)
  }
}

# PHATE -------------------------------------------------------------------

RunPHATE <- function(seu, dims, thread = 1, verbose = F){
  # phate.import <- import("phate")
  # phate.args <- list(data = seu@reductions$pca@cell.embeddings[, paste0("PC_", dims)],
  #                    ndim = 2, knn = 5L, decay = 40, n.landmark = 2000,
  #                    gamma = 1, t = "auto", mds.solver = "sgd", knn.dist.method = "euclidean",
  #                    knn.max = NULL, init = NULL, mds.method = "metric", mds.dist.method = "euclidean",
  #                    t.max = 100, npca = 100, plot.optimal.t = FALSE, verbose = 1,
  #                    n.jobs = 1, seed = NULL, potential.method = NULL, k = NULL, a = NULL)
  # phate <- do.call(what = phate.import$PHATE, args = phate.args)
  #
  # b <- phate.import$PHATE(n_components = 2L, knn = 5L,
  #                         decay = 40, t = "auto", n_landmark = 2000L, mds_solver = "sgd",
  #                         gamma = 1, n_pca = 100L, mds = "metric", mds_dist = "euclidean",
  #                         knn_dist = "euclidean", knn_max = NULL, n_jobs = as.integer(thread),
  #                         random_state = NULL, verbose = F)

  embedding <- suppressWarnings(phateR::phate(data = seu@reductions$pca@cell.embeddings[, paste0("PC_", dims)],
                                              n.jobs = thread)$embedding)
  colnames(embedding) <- c("phate_1", "phate_2")
  rownames(embedding) <- colnames(seu)
  seu@reductions$PHATE <- CreateDimReducObject(
    embeddings = embedding,
    assay = "RNA",
    global = T,
    key = "phate_"
  )
  return(seu)
}
