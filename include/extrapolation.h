#pragma once
#include <cstdlib>
#include <iostream>
#include <cmath>


// For extrapolation
namespace EP {
  // Sequence generator:

  /* function: generate Romberg (even) sequence {h, h/2, h/4, h/8 ...}
     argument: step: sequence array for storing the result
               itermax: array size
  */
  void seq_Romberg (int step[], const std::size_t itermax) {
    if (itermax>0) {
      step[0] = 1;
      for (std::size_t i=1;i<itermax;i++) {
        step[i] = 2*step[i-1];
      }
    }
  }

  /* function: generate Bulirsch & Stoer sequence {h, h/2, h/3, h/4, h/6, h/8 ...}
     argument: step: sequence array for storing the result
               itermax: array size
   */
  void seq_BS (int step[], const std::size_t itermax) {
    if (itermax>0) {
      step[0] = 1;
      std::size_t stepeven = 2; 
      std::size_t stepodd = 3;
      for (std::size_t i=1;i<itermax;i++) {
        if (i%2) {
          step[i] = stepeven;
          stepeven = stepeven*2;
        }
        else {
          step[i] = stepodd;
          stepodd = stepodd*2;
        }
      }      
    }
  }      

  /* function: generate E. Hairer (4k) sequences {h, h/2, h/6, h/10, h/14 ...}
     argument: step: sequence array for storing the result
               itermax: array size
   */
  void seq_Hairer (int step[], const std::size_t itermax) {
    if (itermax>0) {
      step[0] = 1;
      step[1] = 2;
      for (std::size_t i=2;i<itermax;i++) {
        step[i] = step[i-1] + 4;
      }
    }
  }

  // Recursive function for extrapolation
  /* Polynomial_recursion_formula
     function: Using Polynomial function: T_ik = T_i,k-1 + (T_i,k-1 - T_i-1,k-1) / [(h_i-k/h_i)^2 -1]
     argument: ti1k1: t_i-1,k-1
               tik1: t_i,k-1
               hr: h_i-k/h_i
     return:   T_ik
  */
  double polynomial_recursive_formula(const double ti1k1, const double tik1, const double hr) {
    if (hr==1.0) {
      std::cerr<<"Error!: h_i-k/h_i should not be 1.0! (Romberg_recursion_formula)";
      abort();
    }
    return tik1 + (tik1 - ti1k1)/(hr*hr - 1);
  }

  /* Rational_recursion formula
     function: Using rational function: T_ik = T_i,k-1 + (T_i,k-1 - T_i-1,k-1) / {(h_i-k/h_i)^2 * [1- (T_i,k-1 - T_i-1,k-1)/(T_i,k-1 - T_i-1,k-2)] -1}
     argument: ti1k2: T_i-1,k-2
               ti1k1: T_i-1,k-1
               tik1:  T_i,k-1
               hr:    h_i-k/h_i
     return:   T_ik
   */
  double rational_recursive_formula(const double ti1k2, const double ti1k1, const double tik1, const double hr) {
    double dt2 = tik1 - ti1k2;
    double dt1 = tik1 - ti1k1;
    if (dt2==0) {
      if (dt1-dt2==0) return ti1k1;
      else return tik1;
    }
    double dt = dt1/(hr*hr * (1 - dt1/dt2) - 1);
    if (std::isinf(dt)) return tik1;
    else return tik1 + dt;
  }

  // extrapolation iteration function
  /* polynomial extrapolation iteration
    function: iterate T_n,{1..n} based on T_n-1,{1...n-1} and T_n,1
    argument: Tn: Array of T_n-1,{1...n-1}, will be updated to T_n,{1..n}
              Tnew: T_n,1, will be updated to T_n
              step: step sequence
              Tsize: data array size
              n: the new iteration step index (count from 0)
   */
  void polynomial_extrapolation (double** Tn, double* Tnew, const int step[], const std::size_t Tsize, const std::size_t n) {
    for (std::size_t j=0; j<n; j++) {
      //templately storage new results to ttemp
      /*
             T_n-1,j [j]
                       -> T_n,j+1 [n]
              T_n,j [Tnew]
      */
      double hr = (double)step[n]/(double)step[n-j-1];
      for (std::size_t k=0; k<Tsize; k++) Tn[n][k] = polynomial_recursive_formula(Tn[j][k],Tnew[k],hr);
   
      //update j order
      /*
              [j] << T_n,j [Tnew]
      */
      memcpy(Tn[j],Tnew,Tsize*sizeof(double));
   
      //shift temp data to index = n.
      /*
              [Tnew] << T_n,j+1 []
      */
      memcpy(Tnew,Tn[n],Tsize*sizeof(double));
    }
  }

  /* rational extrapolation iteration
    function: iterate T_n,{1..n} based on T_n-1,{1...n-1} and T_n,1
    argument: Tn: Array of T_n-1,{1...n-1}, will be updated to T_n,{1..n}
              Tnew: T_n,1, will be updated to T_n
              Tsize: data array size
              n: the new iteration step index (count from 0)
  */
  void rational_extrapolation (double** Tn, double* Tnew, const int step[], const std::size_t Tsize, const std::size_t n) {
    // additional template storage
    double T1[Tsize];
        
    for (std::size_t j=0; j<n; j++) {
      //templately storage new results to ttemp
      /*
                           T_n-1,j [j]
             T_n-1,j-1 [j-1]          -> T_n,j+1 [n]
                           T_n,j [Tnew]
      */
      double hr = (double)step[n]/(double)step[n-j-1];
      if (j==0) for (std::size_t k=0; k<Tsize; k++) Tn[n][k] = rational_recursive_formula(0,Tn[j][k],Tnew[k],hr);
      else for (std::size_t k=0; k<Tsize; k++) Tn[n][k] = rational_recursive_formula(Tn[j-1][k],Tn[j][k],Tnew[k],hr);

      // update j-1 order
      if (j>0) memcpy(Tn[j-1],T1,Tsize*sizeof(double));

      //storage previous extrapolated data in template position t1
      // [t1] << T_n,j [Tnew]
      memcpy(T1,Tnew,Tsize*sizeof(double));
            
      //shift Tn data to Tnew
      // [Tnew] << T_n,j+1 [n]
      memcpy(Tnew,Tn[n],Tsize*sizeof(double));
    }
  }

  /* error estimation
     function: error calculation based on T_n,n and T_n,n-1 (get maximum error)
     argument: Tn: Tn array
               Tsize: data array size
               n: iteration step index (count from 0)
     return:   maximum error
  */
  double extrapolation_error (double** Tn, const std::size_t Tsize, const std::size_t n) {
    double ermax=0;
    for (std::size_t k=0; k<Tsize; k++) {
      double dik = Tn[n][k];
      double di1k= Tn[n-1][k];
      ermax = std::max(ermax, 2*(dik-di1k)/std::sqrt(dik*dik+di1k*di1k));
    }
    return ermax;
  }

  /* next step optimized factor estimation (based on n)
     function: calculate the modification factor for next extrapolation intergration step (assume next maximum iteration step number is n). Hnew ~ H * (exp/err)**1/(2n+3) 
     argument: err: error of current extrapolation from T_n,n-1 to T_n,n
               exp: expected error
               n: current iteration step index (count from 0)
     return:   optimized factor (H = H*factor)
  */
  double H_opt_factor (const double err, const double exp, const std::size_t n) {
    if (err>0) return std::pow(exp/err, 1/double(2*n+3));
    else return 1.0;
  }
  
}
