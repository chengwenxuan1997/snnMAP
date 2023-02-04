#include <iostream>
#include <queue>
#include <vector>
#include <fstream>
#include <limits>
#include <functional>
// #include <Rcpp.h>
#include <RcppArmadillo.h>
#include <string>

// [[Rcpp::depends(RcppArmadillo)]]


using namespace std;
using namespace Rcpp;

Rcpp::NumericVector CharToNum(CharacterVector Char){
  std::vector<std::string> vstrings(Char.size());
  Rcpp::NumericVector Num(Char.size());
  for (int i=0;i<Char.size();i++){
    vstrings[i]=Char(i);
    Num(i)=std::stoi(vstrings[i]);
  }
  return Num;
}

// Rcpp::CharacterVector SetNames(NumericVector Num){
//   CharacterVector Name(Num.size());
//   for (int i=0;i<Num.size();i++){
//    Name[i]=i+16;
//   }
//   return Name;
// }

Rcpp::NumericVector order_cpp(Rcpp::NumericVector invec){
  int leng = invec.size();
  NumericVector y = clone(invec);
  for(int i=0; i<leng; ++i){
    y[sum(invec<invec[i])] = i+1;
  }
  return(y);
}

void test2(const NumericVector& x) {
  arma::vec yvec(as<arma::vec>(x));
}

Rcpp::NumericVector Rank(arma::vec x) {
  arma::vec test2(x);
  return(Rcpp::as<Rcpp::NumericVector>(wrap(arma::sort_index(arma::sort_index(x))+1)) );
}




// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::export]]

Rcpp::NumericVector DEMaP(std::vector<int> gfrom,
                          std::vector<int> gto,
                          std::vector<double> gw,
                          int NbNodes,
                          std::vector<int> dep, 
                          std::vector<int> arr,
                          NumericMatrix embedding){
  
  //std::vector<std::vector<int> > result;
  Rcpp::NumericVector result(dep.size());

  struct comp{
    bool operator()(const std::pair<int, double> &a, const std::pair<int, double> &b){
      return a.second > b.second;
    }
  };
  
  //Graph
  
  int NbEdges=gfrom.size();
  
  std::vector<std::vector<std::pair<int, double> > > G(NbNodes);                                    
  
  for (int i = 0; i < NbEdges; ++i) {
    G[gfrom[i]].push_back(std::make_pair(gto[i], gw[i]));
  }
  
  //Boucle sur chaque trajet
  std::vector<double> Distances(NbNodes, std::numeric_limits<double>::max());
  for (int j=0; j<dep.size();j++){
    if (j % 256){
      Rcpp::checkUserInterrupt ();
    }
    
    int StartNode=dep[j];
    
    Distances[StartNode] = 0.0;                                                     
    
    priority_queue<std::pair<int, double>, vector<std::pair<int, double> >, comp > Q;
    Q.push(std::make_pair(StartNode, 0.0));                                             
    
    while (!Q.empty()) {                                                          
      int v = Q.top().first;                                                      
      double w = Q.top().second;                                                   
      Q.pop();
      
      if (w <= Distances[v]) {                                                    
        
        for (int i=0; i< G[v].size(); i++) {
          std::pair<int,double> j = G[v][i];                                                  
          int v2 = j.first;                                                      
          double w2 = j.second;
          
          if (Distances[v] + w2 < Distances[v2]) {                                
            Distances[v2] = Distances[v] + w2;                                    
            Q.push(make_pair(v2, Distances[v2]));
          }
        }
      }
    }
    
    Rcpp::NumericVector geo_dist(arr.size());
    for (int i=0;i<arr.size();i++){
      if (Distances[arr[i]]==std::numeric_limits<double>::max()){
        geo_dist[i] = Rcpp::NumericVector::get_na();
      }
      else {
        geo_dist[i] = Distances[arr[i]];
      }
    }
    // std::sort(geo_dist.begin(), geo_dist.end());
    // NumericVector rank = CharToNum(geo_dist.names());
    NumericVector geo_rank = Rank(geo_dist);

    Rcpp::NumericVector eud_dist(arr.size());
    for (int i=0;i<arr.size();i++){
      eud_dist[i] = pow(embedding(i, 0)-embedding(j, 0), 2) + pow(embedding(i, 1)-embedding(j, 1), 2);
    }
    NumericVector eud_rank = Rank(eud_dist);

    double delta = 0;
    for (int i=0;i<arr.size();i++){
      delta += pow(eud_rank[i]-geo_rank[i], 2);
    }
    // result[j] = 1;
    result[j] = 1-6*delta/arr.size()/(pow(arr.size(), 2)-1);
    // result.row(j) = eud_rank;

    //Reinitialize
    std::fill(Distances.begin(),Distances.end(),std::numeric_limits<double>::max());
  }
  
  return result;
  
}


