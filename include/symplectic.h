#pragma once
#include <cstdlib>
#include <iostream>
#include <cmath>

//! For High order Symplectic integrator constructor (Yoshida 1990)
namespace SYM {
    //! Recursive Symplectic cofficients generator 
    /*!
      With Halmitionian Splitting \f$ H = T(p) + U(q) \f$
      The corresponding Symplectic mapping with drift (D) for q and Kick (K) for p is
      \f$ \dot z = \prod_{i=0}^k (exp[c_{k} t D] exp[c_{k+1} t K] exp[c_{k+2} t D]) z \f$
      where \f$ z=(p,q), k = 3^(n-1) \f$
      
      @param[in] coff_array: array to store cofficients \f$ c_k \f$ for steps 
      @param[in] n: half of the order of integrator (order is 2n)
    */
    void recursive_symplectic_cofficients_split(double* coff_array, const int n) {
        if (n==1) {
            coff_array[0] = 0.5;
            coff_array[1] = 1.0;
            coff_array[2] = 0.5;
        }
        else if (n>1) {
            // order 2(n-1) 
            recursive_symplectic_cofficients_split(coff_array, n-1);

            // last order array size
            const int n_size=std::pow(3, n-1);

            // backup last order cofficients
            double coff_backup[n_size];
            for (int i=0; i<n_size; i++) coff_backup[i] = coff_array[i];

            // 2^(1/(2n-1))
            double ap = std::pow(2.0, 1.0/(2*n-1));

            double w[3];
            w[0] = 1.0/(2.0-ap);
            w[1] = - ap*w[0];
            w[2] = w[0];

            for (int i=0; i<3; i++) {
                for (int j=0; j<n_size; j++) {
                    coff_array[i*n_size+j] = w[i]*coff_backup[j];
                }
            }
        }
        else {
            std::cerr<<"Error: integrator order should be positive, given value "<<n<<std::endl;
            abort();
        }
    }

    //! Symplectic cofficients generator
    /*!
      With Halmitionian Splitting \f$ H = T(p) + U(q) \f$
      The corresponding Symplectic mapping with drift (D) for q and Kick (K) for p is
      \f$ \dot z = \prod_{i=0}^k (exp[c_{k} t D] exp[d_{k} t K]) z \f$
      where \f$ z=(p,q), k = 3^(n-1)+1 \f$
      
      call #recursive_symplectic_cofficients_split first, then gether the \f$ c_k, d_k \f$ are gethered in this function
      @params[in] coff_cd: two dimensional cofficients array [1:k][c_k,d_k]
      @params[in] n: half of the order of integrator
      @params[in] n_size: coff_cd array size (for safety check), should be >=k
     */
    void symplectic_cofficients(double coff_cd[][2], const int n, const int c_size) {
        const int n_size = std::pow(3,n);
        double coff_array[n_size];

        recursive_symplectic_cofficients_split(coff_array, n);

        // DKD group number
        const int n_group= n_size/3;

        // coff_cd array size
        const int k = n_group+1;
        
        if(k>c_size) {
            std::cerr<<"Error: coff_cd array size "<<c_size<<" not big enough for order "<<n<<", required size "<<k<<std::endl;
            abort();
        }

        coff_cd[0][0] = coff_array[0];
        for (int i=0; i<n_group-1; i++) {
            coff_cd[i  ][1] = coff_array[3*i+1];
            coff_cd[i+1][0] = coff_array[3*i+2]+coff_array[3*i+3];
        }
        coff_cd[k-2][1] = coff_array[n_size-2];
        coff_cd[k-1][0] = coff_array[n_size-1];
        coff_cd[k-1][1] = 0;
    }
}
