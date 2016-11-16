#pragma once
#include <cstdlib>
#include <iostream>
#include <cmath>
#include <cstring>
#ifdef DEBUG
#include <iomanip>
#endif

//! For extrapolation related functions
namespace EP {
  // Sequence generator:

  //! Generate Harmonic sequence {h, h/2, h/3, h/4 ...}
  /*! @param[out]   step: sequence array for storing the division steps {1, 2, 3, 4 ...}
      @param[in] itermax: array size
   */
  void seq_Harmonic (int step[], const std::size_t itermax) {
    if (itermax>0) {
      step[0] = 1;
      for (std::size_t i=1;i<itermax;i++) {
        step[i] = step[i-1]+1;
      }
    }
  }
  
  //! Generate Romberg (even) sequence {h, h/2, h/4, h/8 ...}
  /*! @param[out]   step: sequence array for storing the division steps {1, 2, 4, 8 ...}
      @param[in] itermax: array size
   */
  void seq_Romberg (int step[], const std::size_t itermax) {
    if (itermax>0) {
      step[0] = 1;
      for (std::size_t i=1;i<itermax;i++) {
        step[i] = 2*step[i-1];
      }
    }
  }

  //! Generate Bulirsch & Stoer sequence {h, h/2, h/3, h/4, h/6, h/8 ...}
  /*! @param[out]   step: sequence array for storing the division steps {1, 2, 3, 4, 6 ...}
      @param[in] itermax: array size
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

  //! Generate E. Hairer (4k) sequences {h/2, h/6, h/10, h/14 ...}
  /*! @param[out]   step: sequence array for storing the division steps {2, 6, 10, 14 ...}
      @param[in] itermax: array size
   */
  void seq_Hairer (int step[], const std::size_t itermax) {
    if (itermax>0) {
      step[0] = 2;
      for (std::size_t i=1;i<itermax;i++) {
        step[i] = step[i-1] + 4;
      }
    }
  }

  // Recursive function for extrapolation
  //! Polynomial_recursion_formula
  /*! Using Polynomial function: \f$ T_{i,k} = T_{i,k-1} + \frac{T_{i,k-1} - T_{i-1,k-1}}{(h_{i-k}/h_i)^2 -1} \f$ \n
     @param [in] ti1k1: \f$ t_{i-1,k-1} \f$
     @param [in] tik1: \f$ t_{i,k-1} \f$
     @param [in] hr: \f$ h_{i-k/h_i} \f$
     \return  \f$ T_{i,k} \f$
  */
  double polynomial_recursive_formula(const double ti1k1, const double tik1, const double hr) {
    if (hr==1.0) {
      std::cerr<<"Error!: h_i-k/h_i should not be 1.0! (Romberg_recursion_formula)";
      abort();
    }
    return tik1 + (tik1 - ti1k1)/(hr*hr - 1);
  }

  //! Rational_recursion formula
  /*! Using rational function: \f$ T_{i,k} = T_{i,k-1} + \frac{ T_{i,k-1} - T_{i-1,k-1} }{ \left( \frac{ h_{i-k} }{ h_i } \right)^2 \left( 1- \frac{ T_{i,k-1} - T_{i-1,k-1} }{ T_{i,k-1} - T_{i-1,k-2} } \right) -1} \f$
     @param [in] ti1k2: \f$ T_{i-1,k-2} \f$ 
     @param [in] ti1k1: \f$ T_{i-1,k-1} \f$ 
     @param [in] tik1:  \f$ T_{i,k-1}   \f$ 
     @param [in] hr:    \f$ h_{i-k/h_i} \f$ 
     \return  \f$ T_{i,k} \f$
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

  // Extrapolation iteration function
  //! Polynomial extrapolation of Tsize number of data together
  /*! Iterate \f$ T_{n,(1..n)} \f$ based on \f$ T_{n-1,(1...n-1)} \f$ and \f$ T_{n,1} \f$
    Notice the Tsize number of data in \a Tn and \a Tnew are independent data that need to be extrapolated individually (each data is an individual \f$ T_{x,x} \f$) 
    @param [in,out] Tn: two dimensional array of \f$ T_{n-1,(1...n-1)} \f$ (size of [n+1][Tsize]; first [] indicate the extrapolation order (1...n), second [] indicate the different data), will be updated to \f$ T_{n,(1..n)} \f$
    @param [in,out] Tnew: one dimensional array of \f$ T_{n,1} \f$ (size of Tsize with different data), will be updated to \f$ T_n \f$
    @param [in] step: step sequence (from sequence generators)
    @param [in] Tsize: data array size
    @param [in] n: the new iteration step index in \a step (count from 0)
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
      std::memcpy(Tn[j],Tnew,Tsize*sizeof(double));
   
      //shift temp data to index = n.
      /*
              [Tnew] << T_n,j+1 []
      */
      std::memcpy(Tnew,Tn[n],Tsize*sizeof(double));
    }
  }

  //! Rational extrapolation of Tsize number of data together
  /*! Iterate \f$ T_{n,(1..n)} \f$ based on \f$ T_{n-1,(1...n-1)} \f$ and \f$ T_{n,1} \f$ \n
    Notice the Tsize number of data in \a Tn and \a Tnew are independent data that need to be extrapolated individually (each data is an individual \f$ T_{x,x} \f$) 
    @param [in,out] Tn: two dimensial array of \f$ T_{n-1,(1...n-1)} \f$ (size of [n+1][Tsize]; first [] indicate the extrapolation order (1...n), second [] indicate the different data) will be updated to \f$ T_{n,(1..n)} \f$
    @param [in,out] Tnew: one dimensional array of \f$ T_{n,1} \f$ (size of Tsize with different data), will be updated to \f$ T_n \f$
    @param [in] step: step sequence (from sequence generators)
    @param [in] Tsize: data array size
    @param [in] n: the new iteration step index in \a step (count from 0)
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
      if (j>0) std::memcpy(Tn[j-1],T1,Tsize*sizeof(double));

      //storage previous extrapolated data in template position t1
      // [t1] << T_n,j [Tnew]
      std::memcpy(T1,Tnew,Tsize*sizeof(double));
            
      //shift Tn data to Tnew
      // [Tnew] << T_n,j+1 [n]
      std::memcpy(Tnew,Tn[n],Tsize*sizeof(double));
    }
  }

  //! Error estimation
  /*! Error calculation based on \f$ T_{n,n} \f$ and \f$ T_{n,n-1} \f$
     Calculate the value of \f$ 2 \frac{T_{n,n} - T_{n,n-1}}{\sqrt{T_{n,n}^2 + T_{n,n-1}^2}} \f$ via looping all Tsize data and select the maximum value as error
     @param [in] Tn: \f$ T_{n,(1..n)} \f$ array (two dimensional array with size [n][Tsize], where first [] indicate the extrapolation order and second [] indicate individual data
     @param [in] Tsize: data array size
     @param [in] n: iteration step index (count from 0)
     \return   maximum error
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

  //! Next step optimized factor estimation (based on the extrapolation order n)
  /*! Calculate the modification factor for next extrapolation intergration step (assume next maximum extrapolation order is n). \n
     (new step size) Hnew ~ (old step size) H * (\a exp/\a err)**1/(2n+3) 
     @param [in] err: error of current extrapolation from \f$ T_{n,n-1} \f$ to \f$ T_{n,n} \f$
     @param [in] exp: expected error
     @param [in] n: current extrapolation order (count from 0)
     return:   optimized factor (H = H*factor)
  */
  double H_opt_factor (const double err, const double exp, const std::size_t n) {
    if (err>0) return std::pow(exp/err, 1/double(2*n+3));
    else return 1.0;
  }

  //! Binomial coefficients generator
  /*! Generate the next binomial sequence (n 1:n) based on (n-1 1:n-1)
     @param [in] bp: n-1 sequence array (size of n-1)
     @param [in] bn: new sequence array (size of n)
     @param [in] n:  new sequence index 
   */
  void binomial_recursive_generator (int *bn, const int* bp, const std::size_t n) {
    if (n>1) {
      bn[0] = 1;
      bn[n-1] = 1;
      // use recursive formula
      for (std::size_t j=1; j<n-1; j++) bn[j] = bp[j-1] + bp[j];
    }
    else bn[0] = 1;
  }

  //! Hermite interpolation coefficients
  /*! Generate Hermite interpolation polynomial coefficients (see example in https://en.wikipedia.org/wiki/Hermite_interpolation)
     @param [out] coff: two dimensional array storing the interpolation coefficients [ndata][\f$\sum_j nlev_j\f$]
     @param [in] x: one dimensional array that store the positions [npoints]
     @param [in] f: two dimensional array that store the f(x) [npoints][ndata]
     @param [in] df: three dimensional array that store the f^(k)(x) [k][npoints][ndata]
     @param [in] ndata: number of different type of data for interpolation
     @param [in] npoints: number of position points.
     @param [in] nlev:  one dimensional array that store the maximum difference level for each position (f^(0)(x) count as 1)
  */
  void Hermite_interpolation_coefficients (double **coff, const double* x, double** f, double ***df, const int ndata, const int npoints, const int* nlev) {
    // get total coefficient number
    int n = 0;
    for (int i=0; i<npoints; i++) n += nlev[i];

    // template data for iteration
    double* dtemp= new double[n];
    // indicator of which point corresponding to in dtemp (-1: need to calculate, >=0; need to set data from points dk)
    int* dk0 = new int[n];
    int* dkn = new int[n];
    
    // loop different type of datt
    for (int i=0; i<ndata; i++) {
      // set initial value
      int joff = 0;
      for (int j=0; j<npoints; j++) {
        for (int k=0; k<nlev[j]; k++) {
          dtemp[k+joff] = f[j][i];
          dk0[k+joff] = j;
          dkn[k+joff] = j;
#ifdef DEBUG
          std::cerr<<std::setprecision(3)<<std::setw(3)<<dk0[k+joff]<<std::setw(3)<<dkn[k+joff]<<std::setw(9)<<dtemp[k+joff]<<"|";
          if (k+joff==n-1) std::cerr<<"\n";
#endif
        }
        joff += nlev[j];
      }

      // first coff (f(x));
      coff[i][0] = dtemp[0];
      // (k!) for cofficients
      int nk=1;
      // j indicate derivate order
      for (int j=1; j<n; j++) {
        nk *=j;
        
        // iteration to get cofficients 
        for (int k=0; k<n-j; k++) {
          // set known derivate value 
          if (dk0[k]==dkn[k+1]) dtemp[k] = df[j-1][dk0[k]][i]/(double)nk;
          // calculate new
          else {
            dtemp[k] = (dtemp[k+1] - dtemp[k])/(x[dkn[k+1]] - x[dk0[k]]);
            dkn[k] = dkn[k+1];
          }
#ifdef DEBUG
          std::cerr<<std::setw(3)<<dk0[k]<<std::setw(3)<<dkn[k]<<std::setw(9)<<dtemp[k]<<"|";
          if (k==n-j-1) std::cerr<<"\n";
#endif
        }

        // get coff
        coff[i][j] = dtemp[0];
      }
    }

    delete[] dtemp;
    delete[] dk0;
    delete[] dkn;
  }

  //! Hermite interpolation polynomial
  /*! Return the interpolation result based on polynomial coefficients generated from EP::Hermite_interpolation_coefficients()
     @param [in] xn: position want to get interpolation value
     @param [out] fxn: one dimensional array that store the interpolation results [ndata]
     @param [in] coff: two dimensional array storing the interpolation coefficients [ndata][\f$\sum_j nlev_j\f$]
     @param [in] x: one dimensional array that store the known positions [npoints].
     @param [in] ndata: number of different type of data for interpolation
     @param [in] npoints: number of position points.
     @param [in] nlev:  one dimensional array that store the maximum difference level for each position (f^(0)(x) count as 1)
  */
  void Hermite_interpolation_polynomial(double xn, double *fxn, double** coff, const double *x, const int ndata, const int npoints, const int* nlev) {
    // offset for shifting dx calculation
    int noff[npoints];
    for (int i=0; i<npoints-1; i++) noff[i] = nlev[i];
    noff[npoints-1] = nlev[npoints-1]-1;

    // polynomial
    for (int i=0; i<ndata; i++) {
      // first cofficient is constant
      fxn[i] = coff[i][0];
      // initial dx
      double dx=1;
      // coff index
      std::size_t ik=1;
      for (int k=0; k<npoints; k++) {
        double dk = xn - x[k];
        for (int j=0; j<noff[k]; j++) {
          dx *= dk;
          fxn[i] += coff[i][ik]*dx;
          ik++;
        }
      }
    }
  }
  
}
