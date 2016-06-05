#include <Rcpp.h>
#include "interaction.h"
#include <cmath>
#include <vector>
#include <string>
#include <sstream>
#include <algorithm>

// [[Rcpp::plugins(cpp11)]]

using namespace std;

double compute_scaling(vector<double> const& attr_data, bool factor, double value){
    double sc = 0;

    if(factor){
        // Inactive fraction
        sc = count_if(attr_data.begin(),attr_data.end(),[value](double a){return a != value;});
        sc /= attr_data.size();
    }
    else{
        for(int i=0;i<attr_data.size();++i){ // average of squared distances
            sc += pow(attr_data[i] - value,2.0);
        }
        sc /= attr_data.size();
    }

    return sc;
}

double InteractionItem::get_value() const{
    return value;
}

string InteractionItem::get_name() const{
  return name;
}

double InteractionItem::get_diff(double other_val) const{
  if(factor){
    if(value == other_val)
      return 0.0;
    else
      return 1.0;
  }
  else{
    return pow(value - other_val,2.0);
  }
}

void InteractionItem::intersect(double data){
    // Compute local error
    double local_error = get_diff(data);

    // Update current error + mark remove if necessary
    current_error += local_error;
    if(current_error > max_error)
        remove = true;
}

string InteractionItem::as_string() const{
  stringstream ss;
  ss << name << "=" << value;
  return ss.str();
}

int InteractionItem::get_idx() const{
    return attr_idx;
}

bool InteractionItem::to_remove() const{
    return remove;
}

InteractionItem::InteractionItem(){
  name = "Default";
  value = 0;
  attr_idx = 0;
  factor = false;
  max_error = 0;
  current_error = 0;
  remove = true;
}

InteractionItem::InteractionItem(string p_name,double p_value,int p_idx,double epsilon,bool p_factor,vector<double> const& attr_data){
    name = p_name;
    value = p_value;
    attr_idx = p_idx;
    factor = p_factor;
    max_error = epsilon*compute_scaling(attr_data,p_factor,p_value);

    current_error = 0;
    remove = false;
}

int Interaction::size() const{
    return items.size();
}

string Interaction::as_string() const{
  stringstream ss;
  
  for(int i=0;i<size();++i){
    ss << items[i].as_string();
    if(i < size()-1) // Not last
      ss << " ";
  }
  
  return ss.str();
}

Rcpp::List Interaction::as_List(Rcpp::List const& datas){
  Rcpp::List inter(size());
  Rcpp::CharacterVector names(size());
  vector<string> info({"attr_idx","value"});
  
  // Interaction
  for(int i=0;i<size();++i){
    names[i] = items[i].get_name();
    
    Rcpp::NumericVector elem(2);
    elem.names() = info;
    elem[0] = items[i].get_idx();
    elem[1] = items[i].get_value();
    
    inter[i] = elem;
  }
  
  inter.names() = names;
  
  // If sims not valid, compute them
  if(!sims_valid)
    compute_sims(datas);
    
  // Add sims
  Rcpp::NumericVector s(sims.size());
  s.names() = datas.names();
  for(int c=0;c<sims.size();++c)
    s[c] = sims[c];
    
  // Result
  vector<string>info2({"interaction","prevalence"});
  Rcpp::List res(2);
  res[0] = inter;
  res[1] = s;
  res.names() = info2;
  
  return res;
}

bool Interaction::check_sims(Rcpp::List const& datas,vector<double> const& theta,bool es){
  if(!es)
    return true;
    
  if(!sims_valid)
    compute_sims(datas);
    
  for(int c=0;c<sims.size();++c){
    if(sims[c] < theta[c])
      return true;
  }
  return false;
}

void Interaction::compute_sims(Rcpp::List const& datas){
  for(int c=0;c<sims.size();++c){ // Iterate over classes
    
    Rcpp::DataFrame class_data = Rcpp::as<Rcpp::DataFrame>(datas[c]);
    int nrows = Rcpp::as<Rcpp::NumericVector>(class_data[0]).size();
    vector<double> diffs(nrows);
    double avg_sim = 0.0;
    
    for(int i=0;i<nrows;++i){ // Interate over instances in class c
      for(int j=0;j<size();++j){ // Iterate over Interaction items
        int attr_idx = items[j].get_idx();
        double instance_val = Rcpp::as<Rcpp::NumericVector>(class_data[attr_idx])[i];
        
        // Update diff
        diffs[i] += items[j].get_diff(instance_val);
      }
      
      avg_sim += exp(-diffs[i]);
    }
    
    avg_sim /= nrows;
    sims[c] = avg_sim;
  }
  
  sims_valid = true;
}

void Interaction::intersect(vector<double> const& instance){
    // Update items
    for(int i=0;i<items.size();++i){
        int idx = items[i].get_idx();
        items[i].intersect(instance[idx]);
    }

    // Remove items
    vector<InteractionItem> new_items(items.size());
    auto it = copy_if(items.begin(),items.end(),new_items.begin(),[](InteractionItem i){return !i.to_remove();});
    new_items.resize(distance(new_items.begin(),it));
    
    // Check sims
    if(items.size() != new_items.size())
      sims_valid = false;
      
    items = new_items;
}

void Interaction::intersect(int row,Rcpp::DataFrame const& data){
    // Update items
    for(int i=0;i<items.size();++i){
        int idx = items[i].get_idx();
        double val = Rcpp::as<Rcpp::NumericVector>(data[idx])[row];
        items[i].intersect(val);
    }

    // Remove items
    vector<InteractionItem> new_items(items.size());
    auto it = copy_if(items.begin(),items.end(),new_items.begin(),[](InteractionItem i){return !i.to_remove();});
    new_items.resize(distance(new_items.begin(),it));
    
    // Check sims
    if(items.size() != new_items.size())
      sims_valid = false;
      
    items = new_items;
}

int Interaction::get_depth() const{
  return depth;
}

void Interaction::set_depth(int i){
  depth = i;
}

Interaction::Interaction(){
  vector<InteractionItem> new_items;
  items = new_items;
  depth = 0;
  sims_valid = false;
  vector<double> new_sims;
  sims = new_sims;
}

Interaction::Interaction(vector<double> const& instance,vector<bool> const& factor,
Rcpp::DataFrame const& class_data,double epsilon_cont,double epsilon_cat,int p_depth,int nb_class){
    vector<InteractionItem> new_items(class_data.size());
    Rcpp::CharacterVector names = class_data.names();

    for(int i=0;i<class_data.size();++i){
        vector<double> attr_data = Rcpp::as<vector<double> >(class_data[i]);
        string name = Rcpp::as<string>(names[i]);
        if(factor[i]){
            InteractionItem item(name,instance[i],i,epsilon_cat,factor[i],attr_data);
            new_items[i] = item;
        }
        else{
            InteractionItem item(name,instance[i],i,epsilon_cont,factor[i],attr_data);
            new_items[i] = item;
        }
    }
    
    items = new_items;
    depth = p_depth;
    sims_valid = false;
    vector<double> new_sims(nb_class);
    sims = new_sims;
}
