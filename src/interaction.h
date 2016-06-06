#ifndef INTERACT_H
#define INTERACT_H

#include <vector>
#include <string>
#include <Rcpp.h>

using namespace std;

double compute_scaling(vector<double> const&, bool, double);

class InteractionItem {
    private:
        string name;
        double value;
        int attr_idx;
        double current_error;
        double max_error;
        bool factor;
        bool remove;

    public:
        void intersect(double);
        int get_idx() const;
        string get_name() const;
        bool to_remove() const;
        double get_value() const;
        double get_diff(double) const;
        string as_string() const;
        InteractionItem();
        InteractionItem(string,double,int,double,bool,vector<double> const&);
};

class Interaction{
    private:
        vector<InteractionItem> items;
        int depth;
        vector<double> sims;
        double radius;
        bool sims_valid;
        void compute_sims(Rcpp::List const&);

    public:
        int size() const;
        void intersect(vector<double> const&);
        void intersect(int,Rcpp::DataFrame const&);
        bool check_sims(Rcpp::List const&,vector<double> const&,bool);
        string as_string() const;
        Rcpp::List as_List(Rcpp::List const&);
        int get_depth() const;
        void set_depth(int);
        Interaction();
        Interaction(vector<double> const&,vector<bool> const&,Rcpp::DataFrame const&,double,double,int,int,double);
};

#endif
