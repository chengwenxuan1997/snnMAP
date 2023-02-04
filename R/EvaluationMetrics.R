# EvaluationMetrics

# if (!grepl("/home/anaconda3/envs/r-4.1.1/bin", Sys.getenv("PATH"))){
#   Sys.setenv("PATH" = paste0(Sys.getenv("PATH"), ":/home/anaconda3/envs/r-4.1.1/bin"))
# }
#
# Rcpp::sourceCpp("/share/genomics/cwx/Research/snnMAP/Codes/EvaluationMetrics.cpp")
# Rcpp::sourceCpp("/share/genomics/cwx/Research/snnMAP/Codes/EvaluationMetricsParallel.cpp")
# Rcpp::sourceCpp("/share/genomics/cwx/Research/snnMAP/Codes/snnMAPhelper.cpp")

# Exponential Attenuation ------------------------------------------------
## turn the fitted formula into function for ggplot

fit2fun <- function(fit.res, var.list = NULL){
  model <- formula(fit.res)
  para.list <- coefficients(fit.res)
  para.list <- c(para.list, var.list)
  curve.formula <- as.character(model)[3]
  for (i in names(para.list)){
    curve.formula <- gsub(i, para.list[i], curve.formula)
  }
  curve.formula <- paste0("function(x) {", curve.formula, "}")
  eval(parse(text = curve.formula))
}


# Calculate DEMaP ---------------------------------------------------------
## denoised embedding manifold preservation

Cal_DEMaP <- function(seu, reduction = 'umap', ncores = min(10, RcppParallel::defaultNumThreads()-2), method = "C"){
  RcppParallel::setThreadOptions(numThreads = ncores)
  cat("prepare the graph \n")
  edges <- Seurat::as.Neighbor(seu@graphs$RNA_nn)@nn.idx
  edges <- as.vector(edges)
  edges <- data.frame("from" = rep(seq(1:ncol(seu)), 20),
                      "to" = edges,
                      "dist" = 1)
  graph <- cppRouting::makegraph(df = edges, directed = F)
  nodes <- unique(c(edges$from, edges$to))

  cat("calculate the distance and correlation coefficient \n")
  if (method == "R"){
    # geo_dist <- as.dist(cppRouting::get_distance_matrix(Graph = graph, from = nodes, to = nodes, allcores = T))
    # embed_dist <- dist(seu@reductions[[reduction]]@cell.embeddings)
    dist <- array(data = c(cppRouting::get_distance_matrix(Graph = graph, from = nodes, to = nodes, allcores = T),
                           as.matrix(dist(seu@reductions[[reduction]]@cell.embeddings))),
                  dim = c(ncol(seu), ncol(seu), 2))
    cat("calculate spearman correlation")
    spearman <- pbapply::pbapply(dist, 1, function(x){
      cor(x[, 1], x[, 2], method = "spearman")
    })
  }else if(method == "C"){
    from = nodes;to = nodes;
    from <- as.character(from)
    to <- as.character(to)
    allnodes <- c(from, to)
    from_id <- graph$dict$id[match(from, graph$dict$ref)]
    to_id <- graph$dict$id[match(to, graph$dict$ref)]
    if (ncores == 1){
      return(DEMaP(gfrom = graph$data[, 1],
                   gto = graph$data[, 2],
                   gw = graph$data[, 3],
                   NbNodes = graph$nbnod,
                   dep = from_id,
                   arr = to_id,
                   embedding = seu@reductions[[reduction]]@cell.embeddings))
    }else{
      return(DEMaP_par(gfrom = graph$data[, 1],
                       gto = graph$data[, 2],
                       gw = graph$data[, 3],
                       NbNodes = graph$nbnod,
                       dep = from_id,
                       arr = to_id,
                       embedding = seu@reductions[[reduction]]@cell.embeddings))
    }
  }
}


# Evaluate the consistency between data and embedding ---------------------


Cal_ARI <- function(seu, reduction = 'umap', method = "hierachy"){
  raw.cluster <- seu$seurat_clusters
  if (method == "hierachy"){
    embedding.cluster <- cutree(tree = hclust(dist(seu@reductions[[reduction]]@cell.embeddings)),
                                k = length(levels(raw.cluster)))
  }else if(method == "kmeans"){
    embedding.cluster <- kmeans(x = seu@reductions[[reduction]]@cell.embeddings,
                                centers = length(levels(raw.cluster)))$cluster
  }

  dict <- as.data.frame(table(raw.cluster, embedding.cluster))
  dict <- dplyr::arrange(dict, -dict$Freq)
  dict <- dict[!duplicated(dict$raw.cluster), ]
  aricode::ARI(c1 = raw.cluster, c2 = embedding.cluster)
}


CalDensityCor <- function(seu, reduction){
  knn = seu@graphs$RNA_nn
  snn = seu@graphs$RNA_snn
  density = rowSums(knn * snn)

  config = uwot::umap.defaults
  arg = uwot:::find_ab_params(config$spread, config$min_dist)
  a = arg["a"]; b = arg["b"]

  radius <- pbapply::pbsapply(1:nrow(snn), function(i){
    neigh <- snn[i, ]
    neigh <- neigh[neigh>0]
    d2 <- as.matrix(dist(seu@reductions[[reduction]]@cell.embeddings[names(neigh), ]))[colnames(seu)[i], names(neigh)]
    q <- 1/(1+a*d2^b)
    sum(q*d2)/sum(q)
  })
  cor(density, radius)
}

CalRadiusCor <- function(seu, reduction, dim.to.use){
  HdRadius <- pbapply::pbsapply(1:nrow(snn), function(i){
    neigh <- snn[i, ]
    neigh <- neigh[neigh>0]
    d2 <- as.matrix(dist(seu@reductions$pca@cell.embeddings[names(neigh), 1:15]))[colnames(seu)[i], names(neigh)]
    q <- 1/(1+a*d2^b)
    sum(q*d2)/sum(q)
  })
  LdRadius <- pbapply::pbsapply(1:nrow(snn), function(i){
    neigh <- snn[i, ]
    neigh <- neigh[neigh>0]
    d2 <- as.matrix(dist(seu@reductions[[reduction]]@cell.embeddings[names(neigh), ]))[colnames(seu)[i], names(neigh)]
    q <- 1/(1+a*d2^b)
    sum(q*d2)/sum(q)
  })
  cor(HdRadius, LdRadius)
}


# Estimate local radius ---------------------------------------------------


CalAttenuation <- function(seu, dims, k.param = 20, complexity = 2, ncores = getOption("mc.cores", 2L), n.iter = 100){
  pcs <- seu@reductions$pca@cell.embeddings[, dims]
  snn <- SeuratObject::as.sparse(Seurat:::FindNeighbors(pcs, k.param = k.param)$snn)
  snn <- as(snn, "dgTMatrix")

  method_args = list("gamma" = 1, "approx_pow" = F)
  # method_args[c("a", "b")] = uwot:::find_ab_params(config$spread, config$min_dist)
  n_epochs = 500; n_components = 2; n_vertices = ncol(seu); negative_sample_rate = 5; grain_size = 1; alpha = 1
  full_opt_args <- uwot:::get_opt_args(opt_args = NULL, alpha)
  epoch_callback = NULL; pcg_rand = T; batch = F; n_sgd_threads = ncores

  message("calculate the graph for ", ncol(seu), " cells ")
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
  message("Fail to estimate the density ", sum(unlist(lapply(fit, is.null))), " cell")

  do.call(rbind, fit)
}

CalRadius <- function(seu, dims, k.param = 20, n_epochs = 500){
  pcs <- seu@reductions$pca@cell.embeddings[, dims]
  knn <- SeuratObject::as.sparse(Seurat::FindNeighbors(seu, dims = dims, k.param = k.param)@graphs$RNA_nn)
  knn <- as(knn, "dgTMatrix")
  knn <- data.frame("from" = knn@i+1 , "to" = knn@j+1)
  knn$dist <- rcpp_distance(mat = pcs, X = knn$from, Y = knn$to)
  knn <- dplyr::arrange(knn, knn$from, knn$dist)
  knn <- split(knn, knn$from)

  smooth.dist <- function(i.dist, iterations = 64, local.connectivity = 1, n.neighbors = Inf,
                          bandwidth = 1, tolerance = 1e-05, min.dist.scale = 0.001){

    i.dist = i.dist[1:min(n.neighbors, nrow(i.dist)), ]
    target = log2(nrow(i.dist)) * bandwidth
    local.int = floor(local.connectivity)
    interpolation = local.connectivity - local.int
    k.dist.mean = mean(i.dist$dist)
    i.dist.cols = seq_len(nrow(i.dist))

    i.nonzero = i.dist$dist[i.dist$dist != 0]
    i.num.nonzero = length(i.nonzero)
    if (i.num.nonzero > local.connectivity){
      rho = i.nonzero[local.int] + interpolation * (i.nonzero[local.int + 1] - i.nonzero[local.int])
    }else{
      rho = max(i.nonzero)
    }

    lo = 0; hi = Inf; mid = 1
    for (n in 2:iterations) {
      val = pmax(i.nonzero - rho, 0)
      val = sum(exp(-val/mid))
      if (abs(val - target) < tolerance) {
        break
      }
      if (val > target) {
        hi = mid
        mid = (lo + hi)/2
      }else {
        lo = mid
        if (!is.finite(hi)) {
          mid = mid * 2
        }else {
          mid = (lo + hi)/2
        }
      }
    }
    result = mid
    result = max(result, min.dist.scale * mean(i.dist$dist))
    cbind(distances = result, nearest = rho)
  }

  knn.smooth <- mclapply(knn, function(knn.dist){
    smooth.dist(i.dist = knn.dist, n.neighbors = k.param)
  })
  knn.smooth <- do.call(rbind, knn.smooth)
  knn.smooth <- list("distances" = knn.smooth[, "distances"],
                   "nearest" = knn.smooth[, "nearest"])

  knn <- do.call(rbind, knn)
  rownames(knn) <- NULL

  bandwidth = 1; mix.ratio = 1
  value <- apply(knn, 1, function(x){
    if(x["dist"] == 0){
      return(0)
    }else{
      return(exp(-(x["dist"] - knn.smooth$nearest[x["from"]])/(knn.smooth$distances[x["from"]] * bandwidth)))
    }
  })

  conn <- cbind(knn, "value" = value)
  conn <- Matrix::sparseMatrix(i = conn$from, j = conn$to, x = conn$value, dims = rep(ncol(seu), 2))
  connP <- conn * Matrix::t(conn)
  graph <- as(mix.ratio * (conn + Matrix::t(conn) - connP) + (1 - mix.ratio) * (connP), "dgTMatrix")
  graph@x[graph@x < max(graph@x)/n_epochs] <- 0
  graph <- as(Matrix::drop0(graph), "dgTMatrix")

  dist <- data.frame("from" = graph@i+1, "to" = graph@j+1)
  dist$dist <- rcpp_distance(mat = pcs, X = dist$from, Y = dist$to)
  dist <- Matrix::sparseMatrix(i = dist$from, j = dist$to, x = dist$dist, dims = rep(ncol(seu), 2))

  rad <- graph*dist*dist
  rad <- rowSums(rad)/rowSums(graph)
}
