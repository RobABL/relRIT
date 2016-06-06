#include <Rcpp.h>
// [[Rcpp::plugins(cpp11)]]

using namespace Rcpp;

// [[Rcpp::export]]
double compute_sim(List const& interaction, DataFrame const& data, LogicalVector const& isFactor, double radius){
  int nrows = as<NumericVector>(data[0]).size();
  int inter_sz = interaction.size();
  NumericVector diffs(nrows);
  double avg_sim = 0.0;

  for(int i = 0;i<nrows;++i){ // iterate over data rows
    for(int j = 0;j<inter_sz;++j){ // iterate over interaction attributes
      List attr_inter = as<List>(interaction[j]);
      int attr_idx = attr_inter["attr_idx"];

      double instance_val = as<NumericVector>(data[attr_idx])[i];
      double inter_val = attr_inter["value"];

      if(isFactor[attr_idx]){ // Categorical attribute
        if(instance_val != inter_val){
          diffs[i]++;
        }
      }
      else{ // Continuous attribute
        diffs[i] += pow(instance_val - inter_val,2.0);
      }
    }
    
    // Sum similarity measures
    avg_sim += std::exp(-diffs[i]/radius);
  }
  
  // Compute mean
  avg_sim /= nrows;
  return avg_sim;
}
