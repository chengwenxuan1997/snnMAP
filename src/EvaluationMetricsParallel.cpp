// [[Rcpp::depends(RcppParallel)]]
// [[Rcpp::plugins(cpp11)]]

#include <iostream>
#include <queue>
#include <vector>
#include <fstream>
#include <limits>
#include <functional>
#include <RcppArmadillo.h>
// #include <Rcpp.h>
#include <string>
#include <RcppParallel.h>

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppParallel)]]
using namespace RcppParallel;
using namespace std;
using namespace Rcpp;

// typedef std::vector<int> stdvec;
// 
// std::vector<double> Rank(arma::vec x) {
//   arma::vec tmp = as<arma::vec>(wrap(arma::sort_index(arma::sort_index(x))+1));
//   // std::vector<int> out=arma::conv_to<stdvec>::from(x);;
//   std::vector<double>out(x.size());
//   return(out);
// }

struct ParDEMaP : public Worker
{
  //input
  const RcppParallel::RVector<int> m_NodeG;
  const RcppParallel::RVector<double> m_WG;
  const RcppParallel::RVector<int> m_IndG;
  RVector<int> m_dep;
  RVector<int> m_arr;
  const int m_nbnodes;
  RcppParallel::RMatrix<double> m_embedding;
  
  //output
  RcppParallel::RVector<double> m_result;
  // RcppParallel::RMatrix<double> m_result;
  
  //constructor
  ParDEMaP(const Rcpp::IntegerVector NodeG,
           const Rcpp::NumericVector WG,
           const Rcpp::IntegerVector IndG,
           Rcpp::IntegerVector dep,
           Rcpp::IntegerVector arr,
           const int nbnodes,
           Rcpp::NumericMatrix embedding,
           Rcpp::NumericVector result) : 
    m_NodeG(NodeG), m_WG(WG), m_IndG(IndG), 
    m_dep(dep), m_arr(arr), m_nbnodes(nbnodes),
    m_embedding(embedding), m_result(result)
  {
    
  }
  
  //overload () operator
  void operator()(std::size_t begin, std::size_t end){
    struct comp{
      
      bool operator()(const std::pair<int, double> &a, const std::pair<int, double> &b){
        return a.second > b.second;
      }
    };
    
    
    std::vector<double> Distances(m_nbnodes, std::numeric_limits<double>::max()); 
    
    for (std::size_t k=begin; k!=end;k++){
      
      int StartNode=m_dep[k];
      
      Distances[StartNode] = 0.0;                                                     
      
      
      
      std::priority_queue<std::pair<int, double>, std::vector<std::pair<int, double> >, comp > Q;
      Q.push(std::make_pair(StartNode, 0.0));                                              
      
      while (!Q.empty()) {                                                        
        int v = Q.top().first;                                                     
        double w = Q.top().second;                                                     
        Q.pop();
        
        if (w <= Distances[v]) {                                                    
          
          for (int i=m_IndG[v]; i< m_IndG[v+1]; i++){
            int v2 = m_NodeG[i];                                                     
            double w2 = m_WG[i];
            
            if (Distances[v] + w2 < Distances[v2]) {                               
              Distances[v2] = Distances[v] + w2;                                   
              
              Q.push(std::make_pair(v2, Distances[v2]));
            }
          }
        }
      }
      
      arma::vec geo_dist(m_nbnodes);
      for (int i=0;i!=m_arr.size();i++){
        if (Distances[m_arr[i]]==std::numeric_limits<float>::max()){
          geo_dist[i]= Rcpp::NumericVector::get_na();
        }
        else {
          geo_dist[i]= Distances[m_arr[i]];
        }
      }
      // std::vector<double> geo_rank = arma::conv_to<std::vector<double>>::from(as<arma::vec>(wrap(arma::sort_index(arma::sort_index(geo_dist))+1)));
      arma::uvec geo_rank = arma::sort_index(arma::sort_index(geo_dist))+1;

      arma::vec eud_dist(m_nbnodes);
      for (int i=0;i<m_nbnodes;i++){
        eud_dist[i] = pow(m_embedding(k, 0)-m_embedding(i, 0), 2) + pow(m_embedding(k, 1)-m_embedding(i, 1), 2);
      }
      arma::uvec eud_rank = arma::sort_index(arma::sort_index(eud_dist))+1;


      double delta = 0;
      for (int i=0;i<m_nbnodes;i++){
        double a = geo_rank(i); double b = eud_rank(i);
        delta += pow(b-a, 2);
      }
      m_result[k] = 1-6*delta/m_nbnodes/(pow(m_nbnodes, 2)-1);
      // m_result[k] = a-b;
      // m_result[k] = eud_rank(0);
      // m_result[k] = k;

      //Reinitialize vectors
      std::fill(Distances.begin(),Distances.end(),std::numeric_limits<double>::max());
      
      
    }
    
  }
  
};


// [[Rcpp::export]]
Rcpp::NumericVector DEMaP_par(Rcpp::IntegerVector gfrom,
                              Rcpp::IntegerVector gto,
                              Rcpp::NumericVector gw,
                              int NbNodes,
                              Rcpp::IntegerVector dep, 
                              Rcpp::IntegerVector arr,
                              Rcpp::NumericMatrix embedding)
  {
  
  
  Rcpp::NumericVector result(NbNodes);
  int NbEdges=gfrom.size();
  
  std::vector<std::vector<std::pair<int, double> > > G(NbNodes);   
  int count=0;
  for (int i = 0; i != NbEdges; ++i) {
    
    G[gfrom[i]].push_back(std::make_pair(gto[i], gw[i]));
    
    count+=1;
  }
  
  
  //Graph vectors
  Rcpp::IntegerVector NodeG(count);
  Rcpp::NumericVector WG(count);
  Rcpp::IntegerVector IndG(NbNodes+1);
  int count2=count;
  count=0;
  for (int i=0; i < G.size();i++){
    IndG[i]=count;
    
    for (int j=0; j < G[i].size();j++){
      NodeG[count]=G[i][j].first;
      WG[count]=G[i][j].second;
      count+=1;
    }
  }
  
  IndG[NbNodes]=count;
  std::vector<std::vector<std::pair<int, double> > > ().swap(G);
  
  
  
  ParDEMaP DEMaPfunc(NodeG,WG,IndG,dep,arr,NbNodes, embedding, result);
  parallelFor(0, dep.length(), DEMaPfunc);
  
  return Rcpp::wrap(result);
  
  
  
}


