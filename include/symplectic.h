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
      @params[in] n: half of the order of integrator, if positive, provide solution one of Yoshida (1990), if negative, provide solution two (only first groups in Table 1/2 of Yoshida 1990), notice only -3 (6th order) and -4 (8th order) works
      @params[in] n_size: coff_cd array size (for safety check), should be >=k
     */
    void symplectic_cofficients(double coff_cd[][2], const int n, const int c_size) {
        if (n>0) {
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
        else if (n==-3) {
            if(c_size<8) {
                std::cerr<<"Error: coff_cd array size "<<c_size<<" not big enough for order "<<n<<", required size "<<8<<std::endl;
                abort();
            }
            coff_cd[0][0] =  0.3922568052387800;
            coff_cd[0][1] =  0.7845136104775600;
            coff_cd[1][0] =  0.5100434119184585;
            coff_cd[1][1] =  0.2355732133593570;
            coff_cd[2][0] = -0.4710533854097566;
            coff_cd[2][1] = -1.1776799841788701;
            coff_cd[3][0] =  0.0687531682525181;
            coff_cd[3][1] =  1.3151863206839063;
            coff_cd[4][0] =  0.0687531682525181;
            coff_cd[4][1] = -1.1776799841788701;
            coff_cd[5][0] = -0.4710533854097566;
            coff_cd[5][1] =  0.2355732133593570;
            coff_cd[6][0] =  0.5100434119184585;
            coff_cd[6][1] =  0.7845136104775600;
            coff_cd[7][0] =  0.3922568052387800;
            coff_cd[7][1] =  0.0000000000000000;
        }
        else if (n==-4) {
            if(c_size<16) {
                std::cerr<<"Error: coff_cd array size "<<c_size<<" not big enough for order "<<n<<", required size "<<16<<std::endl;
                abort();
            }
            /*Solution A
            coff_cd[0][0] =  0.5212131043499550;
            coff_cd[0][1] =  1.0424262086999101;
            coff_cd[1][0] =  1.4313162592035251;
            coff_cd[1][1] =  1.8202063097071399;
            coff_cd[2][0] =  0.9889731189153784;
            coff_cd[2][1] =  0.1577399281236170;
            coff_cd[3][0] =  1.2988836271454836;
            coff_cd[3][1] =  2.4400273261673502;
            coff_cd[4][0] =  1.2164287159851346;
            coff_cd[4][1] = -0.0071698941970812;
            coff_cd[5][0] = -1.2270808589511606;
            coff_cd[5][1] = -2.4469918237052402;
            coff_cd[6][0] = -2.0314077826031052;
            coff_cd[6][1] = -1.6158237415009700;
            coff_cd[7][0] = -1.6983261840452109;
            coff_cd[7][1] = -1.7808286265894515;
            coff_cd[8][0] = -1.6983261840452109;
            coff_cd[8][1] = -1.6158237415009700;
            coff_cd[9][0] = -2.0314077826031052;
            coff_cd[9][1] = -2.4469918237052402;
            coff_cd[10][0] = -1.2270808589511606;
            coff_cd[10][1] = -0.0071698941970812;
            coff_cd[11][0] =  1.2164287159851346;
            coff_cd[11][1] =  2.4400273261673502;
            coff_cd[12][0] =  1.2988836271454836;
            coff_cd[12][1] =  0.1577399281236170;
            coff_cd[13][0] =  0.9889731189153784;
            coff_cd[13][1] =  1.8202063097071399;
            coff_cd[14][0] =  1.4313162592035251;
            coff_cd[14][1] =  1.0424262086999101;
            coff_cd[15][0] =  0.5212131043499550;
            coff_cd[15][1] =  0.0000000000000000;
            */

            /*Solution B*/
            coff_cd[0][0] =  0.4574221231148700;
            coff_cd[0][1] =  0.9148442462297400;
            coff_cd[1][0] =  0.5842687913979845;
            coff_cd[1][1] =  0.2536933365662290;
            coff_cd[2][0] = -0.5955794501471254;
            coff_cd[2][1] = -1.4448522368604799;
            coff_cd[3][0] = -0.8015464361143615;
            coff_cd[3][1] = -0.1582406353682430;
            coff_cd[4][0] =  0.8899492511272584;
            coff_cd[4][1] =  1.9381391376227599;
            coff_cd[5][0] = -0.0112355476763650;
            coff_cd[5][1] = -1.9606102329754900;
            coff_cd[6][0] = -0.9289051917917525;
            coff_cd[6][1] =  0.1027998493919850;
            coff_cd[7][0] =  0.9056264600894914;
            coff_cd[7][1] =  1.7084530707869978;
            coff_cd[8][0] =  0.9056264600894914;
            coff_cd[8][1] =  0.1027998493919850;
            coff_cd[9][0] = -0.9289051917917525;
            coff_cd[9][1] = -1.9606102329754900;
            coff_cd[10][0] = -0.0112355476763650;
            coff_cd[10][1] =  1.9381391376227599;
            coff_cd[11][0] =  0.8899492511272584;
            coff_cd[11][1] = -0.1582406353682430;
            coff_cd[12][0] = -0.8015464361143615;
            coff_cd[12][1] = -1.4448522368604799;
            coff_cd[13][0] = -0.5955794501471254;
            coff_cd[13][1] =  0.2536933365662290;
            coff_cd[14][0] =  0.5842687913979845;
            coff_cd[14][1] =  0.9148442462297400;
            coff_cd[15][0] =  0.4574221231148700;
            coff_cd[15][1] =  0.0000000000000000;
        }
    }
}
