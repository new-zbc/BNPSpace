#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
#include <Rmath.h>
#include <iostream>

using namespace arma;
using namespace Rcpp;
using namespace std;

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends("RcppArmadillo")]]


vec rdirichlet(const vec alpha){
  
  int n = alpha.n_elem;
  vec prob = vec(n, fill::zeros);
  for(int k = 0; k < n; k++){
    prob(k) = as_scalar(randg(1, distr_param(alpha(k), 1.0)));
  }
  
  prob = prob / sum(prob);
  return prob;
}

// calculate the density of Poisson distribution
double dpois(int x, double lambda, bool log_value = true){
  if(lambda <= 0){
    throw std::logic_error("lambda should be positive");
  }
  double r;
  if(x == 0){
    r = - lambda;
  }else{
    vec s1 = arma::linspace(1, x, x);
    r = -sum(log(s1)) + x * log(lambda) - lambda;
  }
  if(log_value){return r;}
  else{return exp(r);}
}


void update_d(mat &D_t, vec &D_0_t, const mat R_t, const uvec gamma_t, const vec group_t, 
              const int K, const mat Y, const vec s, const double a, const double b){
  int J = D_0_t.n_elem;
  
  for(int j = 0; j < J; j++){
    vec y_j = Y.col(j);
    vec R_t_j = R_t.col(j);
    double sample_new = 0;
    
    // cout << "j = " << j << "\n";
    
    if(gamma_t(j) == 0){
      uvec index = find(R_t_j == 0);
      double term1 = a + sum(y_j(index));
      double term2 = b + sum(s(index));
      sample_new = as_scalar(randg(1, distr_param(term1, 1.0/term2)));
      D_0_t(j) = sample_new;
    }
    else{
      for(int k = 0; k < K; k++){
        uvec index = find((R_t_j == 0) && (group_t == k));
        int num = index.n_elem;
        
        // cout << "k = " << k << "\n";
        
        if(num > 0){
          double term1 = a + sum(y_j(index));
          double term2 = b + sum(s(index));
          sample_new = as_scalar(randg(1, distr_param(term1, 1.0/term2)));
        }else{
          sample_new = as_scalar(randg(1, distr_param(a, 1.0/b)));
        }
        D_t(j, k) = sample_new;
      }
    }
  }
}


void update_gamma(uvec &gamma_t, const vec group_t, const mat R_t, const vec s, const int K, 
                  const mat Y, const double aw, const double bw, const double a, const double b){
  int J = gamma_t.n_elem;
  
  uvec index = randi<uvec>(20, distr_param(0, J-1));
  
  
  for(int i = 0; i < 20; i++){
    int j = index(i);
    vec R_t_j = R_t.col(j);
    vec y_j = Y.col(j);
    
    double term1 = 0;
    for(int k = 0; k < K; k++){
      uvec index = find((R_t_j == 0) && (group_t == k));
      int num = index.n_elem;
      
      if(num > 0){
        double a_k = a + sum(y_j(index));
        double b_k = b + sum(s(index));
        term1 += lgamma(a_k) + a * log(b) - lgamma(a) - a_k * log(b_k) + lgamma(aw +1) + lgamma(bw);
      }  
    }
    double term2 = lgamma(a + sum((1-R_t_j) % y_j)) - (a + sum((1-R_t_j) % y_j)) * log(b + sum((1-R_t_j) % s)) + 
      a * log(b) - lgamma(a) + lgamma(aw) + lgamma(bw + 1);
    
    if(term1 - term2 > log(randu())){gamma_t(j) = 1;}
    else{gamma_t(j) = 0;}
  }
}


void update_omega(double &omega_t, const uvec gamma_t, const double aw, const double bw){
  
  int p = gamma_t.n_elem;
  
  double temp = sum(gamma_t);
  
  vec alpha = zeros<vec>(2);
  alpha(0) = aw + temp;
  alpha(1) = bw + p - temp;
  
  omega_t = rdirichlet(alpha)(0);
  
}


void update_R(mat &R_t, const uvec gamma_t, const vec group_t, const int K, 
              const mat D_t, const vec D_0_t, const vec s, const mat Y, const double pi_r){
  int J = gamma_t.n_elem;
  int N = group_t.n_elem;
  
  for(int i = 0; i < N; i++){
    // cout << "i = " << i << "\n";
    
    for(int j = 0; j < J; j++){
      if(Y(i, j) == 0){
        if(gamma_t(j) == 0){
          double log_prob = log(pi_r) - log(pi_r + (1-pi_r)*exp(-s(i) * D_0_t(j)));
          if(log_prob > log(randu())){R_t(i, j) = 1;}
          else{R_t(i, j) = 0;}
        }
        else{
          int k = group_t(i);
          double log_prob = log(pi_r) - log(pi_r + (1-pi_r)*exp(-s(i) * D_t(j, k)));
          if(log_prob > log(randu())){R_t(i, j) = 1;}
          else{R_t(i, j) = 0;}
        }
      }
    }
  }
}


void update_pi(double &pi_t, const mat R_t, const double ar, const double br){
  
  int N = R_t.n_rows;
  int p = R_t.n_cols;
  
  double temp = accu(R_t);
  
  vec alpha = zeros<vec>(2);
  alpha(0) = ar + temp;
  alpha(1) = br + N * p - temp;
  
  pi_t = rdirichlet(alpha)(0);
  
}


//function to calculate MRF energy function
double MRF_energy(const int &k, const int &l, const vec &group_t, const umat &G,
                  const double f, const bool log_value = true){
  
  uvec G_l = G.row(l).t();
  uvec neighbor = G_l.elem(find(G_l > 0)) - 1;
  
  uvec neighbor_id = find(group_t(neighbor) == k);
  double result = f * neighbor_id.n_elem;
  
  if(log_value){return result;}
  else{ return exp(result); }
}


// function to calculate zero-inflated Poisson model
double dZIpoi(const vec y, const vec r, const uvec gamma, const vec d_k, const vec d_0, const double s){
  int J = y.n_elem;
  double result = 0;
  for(int j = 0; j < J; j++){
    if(r(j) == 0){
      if(gamma(j) == 1){result += dpois(y(j), s*d_k(j));}
      else{result += dpois(y(j), s*d_0(j));}
    }
  }
  return result;
}


//function to calculate probability of existing clusters
double prob_existing(const int k, const int l, const vec group_t, const uvec gamma_t, const vec d_k, const vec d_0, 
                     const vec r, const vec y, const double s, const umat G, const double f){
  double result = 0;
  
  result +=  MRF_energy(k, l, group_t, G, f);
  
  result += dZIpoi(y, r, gamma_t, d_k, d_0, s);
  
  return result;
}


void update_group(vec &group_t,  const mat &D_t, const uvec gamma_t, const vec D_0, const mat R_t, const vec s,
                  const mat &Y, const umat &G, const double &f, const double a, const double b, const double &temperature){
  
  int N = Y.n_rows;
  vec clusters_name = unique(group_t);
  int K = clusters_name.n_elem;
  
  uvec S = randperm(N);
  
  
  for(int id = 0; id < N; id++){
    int i = S(id);
    int label_old = group_t(i);
    int label_new;
    uvec group_i = find(group_t == label_old);
    int c_size = group_i.n_elem;
    
    // cout << "i = " << i << "\n";
    
    vec prob_i = zeros<vec>(K);
    for(int k = 0; k < K; k++){
      prob_i(k) = prob_existing(k, i, group_t, gamma_t, D_t.col(k), D_0, R_t.row(i).t(), Y.row(i).t(), s(i), G, f);
    }
      
    prob_i = (prob_i - max(prob_i))/temperature;
    prob_i = exp(prob_i) / sum(exp(prob_i));
      
    label_new = as_scalar(Rcpp::RcppArmadillo::sample(regspace(0, K-1), 1, false, prob_i));
      
    // cout << "label_new = " << label_new << "\n";
    if(c_size > 1){group_t(i) = label_new;}
  }
}


// [[Rcpp::export]]
Rcpp::List run_fix_clusters(vec &group_t,  mat &D_t, uvec &gamma_t, vec &D_0_t, mat &R_t, double &pi_t,
                   const double f, const double a, const double b, const double aw, 
                   const double bw, const double ar, const double br, 
                   const mat Y, const umat G, const vec s, const int max_iters, const int seed=12569){
  int N = Y.n_rows;
  int p = Y.n_cols;
  vec clusters = unique(group_t);
  int K = clusters.n_elem;
  
  //store the parameters
  mat group_iter = zeros<mat>(N, max_iters);
  cube D_iter = zeros<cube>(p, K, max_iters);
  mat D_0_iter = zeros<mat>(p, max_iters);
  umat gamma_iter = zeros<umat>(p, max_iters);
  vec pi_iter = zeros<vec>(max_iters);
  
  //store the initialized pars
  
  group_iter.col(0) = group_t;
  D_iter.slice(0) = D_t;
  D_0_iter.col(0) = D_0_t;
  gamma_iter.col(0) = gamma_t;
  pi_iter(0) = pi_t;
  
  for(int t = 1; t < max_iters; t++){
    
    double temperature = pow(0.999, t);
    
    update_gamma(gamma_t, group_t, R_t, s, K, Y, aw, bw, a, b);
    
    update_d(D_t, D_0_t, R_t, gamma_t, group_t, K, Y, s, a, b); 
    
    update_R(R_t, gamma_t, group_t, K, D_t, D_0_t, s, Y, pi_t);
    
    update_pi(pi_t, R_t, ar, br);
    
    update_group(group_t, D_t, gamma_t, D_0_t, R_t, s, Y, G, f, a, b, temperature);
    
    //store the pars
    group_iter.col(t) = group_t;
    D_iter.slice(t).head_cols(K) = D_t;
    D_0_iter.col(t) = D_0_t;
    gamma_iter.col(t) = gamma_t;
    pi_iter(t) = pi_t;
  }
  
  return List::create(Named("group_iter") = group_iter, 
                      Named("D_iter") = D_iter,
                      Named("D_0_iter") = D_0_iter,
                      Named("gamma_iter") = gamma_iter,
                      Named("pi_iter") = pi_iter);
  
}



