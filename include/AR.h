#pragma once

#include <iostream>
#include <cstring>
#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <limits>
#ifdef USE_OMP
#include <omp.h>
#endif

#ifdef DEBUG
#include <iomanip>
#define WIDTH 10
#endif

#include "extrapolation.h"

//! Algorithmic regularization chain (ARC) namespace
/*!
  All major ARC classes and related acceleration functions (typedef) are defined
 */
namespace ARC {

#ifdef TIME_PROFILE
#include <sys/time.h>
static double get_wtime(){
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return tv.tv_sec + 1.e-6 * tv.tv_usec;
}

//! A structure storing time profile 
class chainprofile{
public:
  double t_apw;    ///< APW calculation
  double t_uplink; ///< update link
  double t_lf;     ///< leap-frog
  double t_ep;     ///< extrapolation
  double t_pext;   ///< perturber force
  double t_dense;  ///< dense output
  double t_init;   ///< initialization

  int* stepcount;  ///< iteration step count in extrapolation
  int  itercount;  ///< total substep in extrapolation

  chainprofile() {reset_tp();}  ///< initialization 

  void initstep(const std::size_t n) {
    if (!stepcount) {
      stepcount=new int[n];
      for (std::size_t i=0; i<n; i++) stepcount[i]=0;
    }
  }

  /*! reset the time*/
  void reset_tp(){
    t_apw=0.0;
    t_uplink=0.0;
    t_lf=0.0;
    t_ep=0.0;
    t_pext=0.0;
    t_dense=0.0;
    t_init=0.0;
    stepcount=NULL;
    itercount=0;
  }

  ~chainprofile() { if (stepcount) delete[] stepcount;}

};

#endif

//declaration
template <class particle, class int_par> class chain;
template <class particle, class int_par> class chainlist;
template <class int_par> class chainpars;


//! The chain parameter controller class
/*!
  This class control the acceleration function and time transformation function (::ARC::chainpars.pair_AW, ::ARC::chainpars.pair_Ap), integration methods and error parameters for \ref chain, and timescale calculator ::ARC::chainpars.pair_T.
  The template depended type int_par is user-defined class for two-body interations used in ::ARC::chainpars.pair_AW, ::ARC::chainpars.pair_Ap and ::ARC::chainpars.pair_T
 */
template <class int_par>
class chainpars{
template <class T, class P> friend class chain;
private:
  // time step integration parameter
  double alpha; ///< logH cofficient
  double beta;  ///< TTL cofficient
  double gamma; ///< constant

  // extrapolation control parameter
  int exp_method;         ///< 1: Polynomial method; others: Rational interpolation method
  int exp_sequence;       ///< 1: Romberg sequence {h, h/2, h/4, h/8 ...}; 2: Bulirsch & Stoer (BS) sequence {h, h/2, h/3, h/4, h/6, h/8 ...}; 3: 4k sequence {h, h/2, h/6, h/10, h/14 ...}; others: Harmonic sequence {h, h/2, h/3, h/4 ...} (defaulted 2. BS sequence)
  std::size_t exp_itermax; ///< maximum times for iteration.

  int* step; ///< substep sequence
  int** bin_index; ///< binomial coefficients

  int auto_step;           ///< if 0: no auto step; if 1: use extrapolation error to estimate next step modification factor; if 2: use min X/(g V) and V/(g A) to obtain next step; if 3: use maximum sequence index to control next step; if 4: use mimimum kepler period of every two neigbor members to estimate the next step modification factor.
  double auto_step_fac_min; ///< minimum reduction factor
  double auto_step_fac_max; ///< maximum reduction factor
  double auto_step_eps;     ///< coefficient
  std::size_t auto_step_iter_min;   ///< mimimum iteration level
  std::size_t auto_step_iter_max;   ///< maximum iteration level

public:
  //! Function pointer type to function for calculation of acceleration (potential) and the component of time transformation function \f$\partial W_{ij}/\partial \mathbf{x}_i\f$ and \f$W_{ij}\f$ (from particle j to particle i).
  /*!          @param[out] Aij: acceleration vector for particle i.
    @param[out] Pij: potential of particle i from particle j (cumulative total potential only count the case of j>i)
    @param[out] dWij: time transformation function partial derivate \f$\partial W_{ij}/\partial \mathbf{x}_i\f$ (component from j to i; used when #beta>0; see ARC::chainpars.setabg()) .
    @param[out] Wij: time transformation function component from j to i (used when #beta>0; notice in \ref ARC::chain the cumulative total W only count the case of j>i).
    @param[in] Xij: relative position (1:3) \f$ \mathbf{x}_j - \mathbf{x}_i \f$
    @param[in] mi: particle i mass.
    @param[in] mj: particle j mass.
    @param[in] pars: User-defined interaction parameter class type data
  */
  typedef void (*pair_AW) (double*, double &, double *, double &, const double*, const double &, const double &, const int_par*);

  //! Function pointer type to function for calculation of acceleration (potential) from particle j to particle i.
  /*!         @param[out]  Aij: acceleration vector of particle i.
    @param[out]  Pij: potential of particle i from j
    @param[in]  pi: position vector i.
    @param[in]  pj: position vector j.
    @param[in]  mi: particle mass i.
    @param[in]  mj: particle mass j.
    @param[in]  pars: User-defined interaction parameter class type data
  */
  typedef void (*pair_Ap) (double *, double &, const double*, const double*, const double&, const double&, const int_par*);

  //! Function pointer type to function for calculation the timescale of two-body motion
  /*!
    @param[in] m1: mass of particle 1
    @param[in] m2: mass of particle 2
    @param[in] X: relative position vector
    @param[in] V: relative velocity vector
    @param[in] pars: User-defined interaction parameter class type data
  */
  typedef double (*pair_T) (const double, const double, const double*, const double*, const int_par*);

  // pair force
  pair_AW pp_AW;  ///< acceleration and dW/dr of two particle
  pair_Ap pp_Ap;  ///< accelaration of two particle
  pair_T  pp_T;   ///< two-body timescale calculator
  
  // time step
  double dtmin; ///< minimum physical time step
  double dterr; ///< physical time error criterion

  // extrapolation control parameter
  double exp_error;        ///< relative error requirement for extrapolation
  bool exp_fix_iter;       ///< flag showing whether the times of iteration is fixed or not

  //! constructor with defaulted parameters
  /*! - Acceleration function use ARC::Newtonian_AW() and ARC::Newtonian_Ap().
      - ARC method use logarithmic Hamiltonian (logH) (#alpha = 1.0, #beta = 0.0 #gamma = 0.0).
      - Phase/energy error limit #exp_error = 1e-10.
      - Minimum physical time step #dtmin = 5.4e-20.
      - Time synchronization error limit #dterr = 1e-10.
      - Maximum extrapolation sequence index (accuracy order/iteration times) #exp_itermax = 20
      - The maximum sequence index (iteration times) is adjusted by error criterion
      - Bulirsch & Stoer sequence {h, h/2, h/3, h/4, h/6...} is used
      - No auto-step
   */
  chainpars(): alpha(1.0), beta(0.0), gamma(0.0), pp_AW(NULL), pp_Ap(NULL), pp_T(NULL) {
    step = NULL;
    bin_index = NULL;
    setEXP(1E-10, 5.4E-20, 1E-10, 20, 2, 2, false);
    setAutoStep(0);
  }

  //! constructor
  /*!
    All parameters can be set except anto-step method (defaulted is zero)
    @param [in] aw: acceleration and dW/dr of two particle calculation function pointer with ::ARC::pair_AW type. (interaction between memebrs)
    @param [in] ap: acceleration calculation function pointer with ::ARC::pair_Ap type. (interaction between member and perturber)
    @param [in] at: two-body timescale calculation function pointer with ::ARC::pair_T type. (this will be used for new step size estimation)
    @param [in] a,b,g: ARC time transformation method coefficients (\f$ dt = ds/[a *(logH) + b * (TTL) + g])\f$. \n
                - a: Logarithmic Hamiltonian (logH) method coefficient #alpha (0.0, 1.0)
                - b: Time-Transformed Leapfrog (TTL) method coefficient #beta (0.0, 1.0)
                - g: Constant coefficient (no time transformation) #gamma
    @param [in] error: Phase/energy error limit (defaulted #exp_error = 1e-10)
    @param [in] dtm: Minimum physical time step (defaulted #dtmin = 5.4e-20)
    @param [in] dte: Time synchronization error limit (defaulted #dterr = 1e-6)
    @param [in] itermax: Maximum extrapolation sequence index (accuracy order/iteration times) (defaulted #exp_itermax = 20)
    @param [in] ext_method: 1: Polynomial interpolation method; others: Rational interpolation method (defaulted: Rational)
    @param [in] ext_sequence: 1: Romberg sequence {h, h/2, h/4, h/8 ...}; 2: Bulirsch & Stoer (BS) sequence {h, h/2, h/3, h/4, h/6, h/8 ...}; 3: 4k sequence {h, h/2, h/6, h/10, h/14 ...}; others: Harmonic sequence {h, h/2, h/3, h/4 ...} (defaulted 2. BS sequence)
    @param [in] ext_iteration_const: true: the maximum sequence index (iteration times) is fixed to itermax; false: adjust the maximum sequence index by error criterion (false)
   */
  chainpars(pair_AW aw, pair_Ap ap, pair_T at, const double a, const double b, const double g, const double error=1E-10, const double dtm=5.4e-20, const double dte=1e-6, const std::size_t itermax=20, const int ext_method=2, const int ext_sequence=2, const bool ext_iteration_const=false) {
    step = NULL;
    bin_index = NULL;
    setabg(a,b,g);
    setEXP(error,dtm,dte,itermax,ext_method,ext_sequence,ext_iteration_const);
    setA(aw,ap,at);
    setAutoStep(0);
  }

  //! destructor
  ~chainpars() {
    if (step!=NULL) delete[] step;
    if (bin_index!=NULL) {
      for (std::size_t i=0;i<exp_itermax;i++) 
        if (bin_index[i]!=NULL)
          delete[] bin_index[i];
      delete[] bin_index;
    }
  }

  //! Set acceleration, potential, time transformation function \f$\partial W/\partial r\f$ and \f$W\f$ calculator
  /*! Set acceleration, potential, time transformation function \f$\partial W/\partial r\f$ and \f$W\f$ for two particles. This is the basic functions used for interaction between chain members and from perturbers.
    @param [in] aw: acceleration, potential and time transformation function \f$\partial W/\partial r\f$, \f$W\f$ of two particles calculation function pointer with ::ARC::pair_AW type. (interaction between memebrs). Notice this function defines \f$\partial W/\partial r\f$ and \f$W\f$, and these time transformation functions are only used when #beta>0 (see setabg()).
    @param [in] ap: acceleration calculation function pointer with ::ARC::pair_Ap type. (interaction between members and perturbers)
    @param [in] at: two-body timescale calculation function pointer with ::ARC::pair_T type. (this will be used for new step size estimation)
  */
  void setA(pair_AW aw, pair_Ap ap, pair_T at=NULL) {
    pp_AW = aw;
    pp_Ap = ap;
    pp_T  = at;
    // safety check
    if (pp_AW==NULL) {
      std::cerr<<"Error: accelaration and time trasnformation calculator function is NULL\n";
      abort();
    }
    if (pp_Ap==NULL) {
      std::cerr<<"Error: accelaration calculator function is NULL\n";
      abort();
    }
  }
  
     
  //! Set time step transformation parameter
  /*!
      Set parameter a,b,g where physical time step \f$ dt = ds/[a *(logH) + b * (TTL) + g]\f$ \n
      @param [in] a: Logarithmic Hamiltonian (logH) method coefficient #alpha (0.0, 1.0)
      @param [in] b: Time-Transformed Leapfrog (TTL) method coefficient #beta (0.0, 1.0)
      @param [in] g: Constant coefficient (no time transformation) #gamma
  */
  void setabg(const double a=1.0, const double b=0.0, const double g=0.0) {
    alpha = a;
    beta = b;
    gamma = g;
    // safety check
    if (alpha==0&&beta==0&&gamma==0) {
      std::cerr<<"Error: alpha, beta and gamma cannot be all zero!\n";
      abort();
    }
  }

  //! Get time step transformation parameter
  /*!
      Get parameter a,b,g where physical time step \f$ dt = ds/[a *(logH) + b * (TTL) + g]\f$ \n
      @param [in] a: Logarithmic Hamiltonian (logH) method coefficient #alpha (0.0, 1.0)
      @param [in] b: Time-Transformed Leapfrog (TTL) method coefficient #beta (0.0, 1.0)
      @param [in] g: Constant coefficient (no time transformation) #gamma
  */
  void getabg(double &a, double &b, double &g) const {
    a = alpha;
    b = beta;
    g = gamma;
  }
  

  //! Set extrapolation parameters
  /*! Set extrapolation related parameters as following:
    @param [in] error: phase/energy relative error requirement for extrapolation (defaulted #exp_error = 1e-10)
    @param [in] dtm: minimum physical time step allown (defaulted #dtmin = 5.4e-20)
    @param [in] dte: Time synchronization error limit (defaulted #dterr = 1e-6)
    @param [in] itermax: Maximum extrapolation sequence index (accuracy order/iteration times) (defaulted #exp_itermax = 20)
    @param [in] method: 1: Polynomial method; others: Rational interpolation method (defaulted Rational)
    @param [in] sequence: 1: Romberg sequence {h, h/2, h/4, h/8 ...}; 2: Bulirsch & Stoer (BS) sequence {h, h/2, h/3, h/4, h/6, h/8 ...}; 3: 4k sequence {h, h/2, h/6, h/10, h/14 ...}; others: Harmonic sequence {h, h/2, h/3, h/4 ...} (defaulted 2. BS sequence)
    @param [in] ext_iteration_const: true: the maximum sequence index (iteration times) is fixed to itermax; false: adjust the maximum sequence index by error criterion (false)
  */
  void setEXP(const double error=1E-10, const double dtm=5.4e-20, const double dte=1e-6, const std::size_t itermax=20, const int method=2, const int sequence=2, const bool ext_iteration_const=false) {
    exp_error = error;
    exp_method = method;
    exp_sequence = sequence;
    exp_fix_iter = ext_iteration_const;
    dterr = dte;
    dtmin = dtm;

    // delete binomial array
    if (bin_index!=NULL) {
      if (step!=NULL) {
        for (std::size_t i=0; i<(std::size_t)step[exp_itermax]; i++)
          if (bin_index[i]!=NULL) delete[] bin_index[i];
        delete[] bin_index;
      }
    }
    
    // reset step array
    if (step!=NULL) delete[] step;
    step = new int[itermax+1];

    // calculate sequences of steps
    // Romberg (even) sequence {h, h/2, h/4, h/8 ...}
    if (sequence==1) EP::seq_Romberg(step,itermax+1);
    // Bulirsch & Stoer sequence {h, h/2, h/3, h/4, h/6, h/8 ...}
    else if (sequence==2) EP::seq_BS(step,itermax+1);
    // E. Hairer (4k) sequences {h, h/2, h/6, h/10, h/14 ...}
    else if (sequence==3) EP::seq_Hairer(step,itermax+1);
    // Harmonic sequences {h, h/2, h/3, h/4 ...}
    else EP::seq_Harmonic(step,itermax+1);

    // calculate binomial coefficients
    bin_index = new int*[step[itermax]];
    for (std::size_t i=0; i<(std::size_t)step[itermax]; i++) {
      bin_index[i] = new int[i+1];
      if (i>0) EP::binomial_recursive_generator(bin_index[i],bin_index[i-1],i+1);
      else EP::binomial_recursive_generator(bin_index[i],NULL,i+1);
    }
    exp_itermax = itermax;
  }

  //! Get interpolation method index
  /*! 
    \return method index: 1: Polynomial method; others: Rational interpolation method (defaulted Rational)
  */
  const int getEXPmethod() const {
    return exp_method;
  }

  //! Get sequence method indices
  /*  
    \return sequence method index: 1: Romberg sequence {h, h/2, h/4, h/8 ...}; 2: Bulirsch & Stoer (BS) sequence {h, h/2, h/3, h/4, h/6, h/8 ...}; 3: 4k sequence {h, h/2, h/6, h/10, h/14 ...}; others: Harmonic sequence {h, h/2, h/3, h/4 ...} (defaulted 2. BS sequence)
  */  
  const int getEXPsequence() const {
    return exp_sequence;
  }

  //! Get extrapolation maximum iteration times (maximum sequence index)
  /*  
    \return maximum iteration times
  */  
  const int getEXPitermax() const {
    return exp_itermax;
  }
  
  //! Determine auto-step parameter
  /*! @param[in] option: auto-step method
    - 0. no auto-step
    - 1. use extrapolation error to estimate next step
    - 2. use min(\f$X/(gV)\f$,\f$V/(gA)\f$) to estimate next step
    - 3: use maximum sequence index to control next step. If index >\a iter_max, modification factor = \a factor_min; if index < \a iter_min, modification factor = \a factor_max;
    - 4: use mimimum user-defined two-body timescale of every two neigbor members to estimate next step
      @param[in] factor_min: if \a option = 1->3, it is the minimum step reduction factor (defaulted: 0.7)
      @param[in] factor_max: if \a option = 1->3, it is the maximum step reduction factor (defaulted: 1.3)
      @param[in] eps: coefficient. If \a option = 2 or 4, coefficient is the multiplying factor for the estimated next step (defaulted: 0.125)
      @param[in] iter_min: if \a option = 3, it is the minimum sequence index (iteration times) criterion during extrapolation intergration (defaulted: 5)
      @param[in] iter_max: if \a option = 3, it is the maximum sequence index (iteration times) criterion during extrapolation intergration (defaulted: 17)
  */
  void setAutoStep(const int option, const double factor_min=0.7, const double factor_max=1.3, const double eps=0.125, const std::size_t iter_min=5, const std::size_t iter_max=17) {
    if (option<0||option>4) {
      std::cerr<<"Error: autostep options should be set from 0 to 4, current value: "<<option<<"!\n";
      abort();
    }
    auto_step = option;
    if (factor_min>1.0||factor_min<0.0) {
      std::cerr<<"Error: step reduction factor is limited to the range [0,1.0], input value is "<<factor_min<<"!\n";
      abort();
    }
    auto_step_fac_min = factor_min;
    if (factor_max<1.0) {
      std::cerr<<"Error: step increasing factor is limited to the range (>1.0), input value is "<<factor_max<<"!\n";
      abort();
    }
    auto_step_fac_max = factor_max;
    if (eps<0.0) {
      std::cerr<<"Error: multiplying coefficient for next estimated step should be positive, input value is "<<eps<<"!\n";
      abort();
    }
    auto_step_eps=eps;
    if (iter_min>exp_itermax) {
      std::cerr<<"Error: minimum iteration times criterion should be less than extrapolation maximum iteration times ("<<exp_itermax<<"), input value is "<<iter_min<<"!\n";
      abort();
    }
    auto_step_iter_min=iter_min;
    if (iter_max>exp_itermax) {
      std::cerr<<"Error: maximum iteration times criterion should be less than extrapolation maximum iteration times ("<<exp_itermax<<"), input value is "<<iter_max<<"!\n";
      abort();
    }
    auto_step_iter_max=iter_max;
  }

  //! Get auto-step parameters
  /*! @param[in] option: auto-step method
    - 0. no auto-step
    - 1. use extrapolation error to estimate next step
    - 2. use min(\f$X/(gV)\f$,\f$V/(gA)\f$) to estimate next step
    - 3: use maximum sequence index to control next step. If index >\a iter_max, modification factor = \a factor_min; if index < \a iter_min, modification factor = \a factor_max;
    - 4: use mimimum user-defined two-body timescale of every two neigbor members to estimate next step
      @param[in] factor_min: if \a option = 1->3, it is the minimum step reduction factor (defaulted: 0.7)
      @param[in] factor_max: if \a option = 1->3, it is the maximum step reduction factor (defaulted: 1.3)
      @param[in] eps: coefficient. If \a option = 2 or 4, coefficient is the multiplying factor for the estimated next step (defaulted: 0.125)
      @param[in] iter_min: if \a option = 3, it is the minimum sequence index (iteration times) criterion during extrapolation intergration (defaulted: 5)
      @param[in] iter_max: if \a option = 3, it is the maximum sequence index (iteration times) criterion during extrapolation intergration (defaulted: 17)
   */
  void getAutoStep(int &option, double &factor_min, double &factor_max, double &eps, std::size_t &iter_min, std::size_t &iter_max) const {
    option = auto_step;
    factor_min = auto_step_fac_min;
    factor_max = auto_step_fac_max;
    eps = auto_step_eps;
    iter_min = auto_step_iter_min;
    iter_max = auto_step_iter_max;
  }

  //! parameters dumping
  /*! Dump parameters into file. Notice dynamic array data will not be dumped.
    The dumping data includes: #alpha, #beta, #gamma, #exp_method, #exp_sequence, #exp_itermax, #exp_error, #exp_fix_iter, #dtmin, #dterr, #auto_step, auto_step parameters
    @param [in] filename: file for storing the data
   */
  void dump(const char* filename) const {
    std::FILE* pout = std::fopen(filename,"w");
    if (pout==NULL) std::cerr<<"Error: filename "<<filename<<" cannot be open!\n";
    else {
      fwrite(&alpha, sizeof(double),1,pout);
      fwrite(&beta,  sizeof(double),1,pout);
      fwrite(&gamma, sizeof(double),1,pout);
      fwrite(&exp_method,  sizeof(int),1,pout);
      fwrite(&exp_sequence,sizeof(int),1,pout);
      fwrite(&exp_itermax, sizeof(std::size_t),1,pout);
      fwrite(&exp_error,   sizeof(double),1,pout);
      fwrite(&exp_fix_iter,sizeof(bool),1,pout);
      fwrite(&dtmin, sizeof(double),1,pout);
      fwrite(&dterr, sizeof(double),1,pout);
      fwrite(&auto_step, sizeof(int),1,pout);
      fwrite(&auto_step_fac_min, sizeof(double),1,pout);
      fwrite(&auto_step_fac_max, sizeof(double),1,pout);
      fwrite(&auto_step_eps,     sizeof(double),1,pout);
      fwrite(&auto_step_iter_min,sizeof(std::size_t),1,pout);
      fwrite(&auto_step_iter_max,sizeof(std::size_t),1,pout);
    }
  }

  //! parameters loading from a dumped file
  /*! Load parameters from a dumped file. All dynamical array will be initialized after reading parameters.
    Notice this function will load all parameters except acceleration function pointers #pp_AW and #pp_Ap. 
    setA() should be called before chain integration.
    The reading data list is shown in dump()
    @param [in] filename: file to read the data
   */
  void load(const char* filename) {
    std::FILE* pin = std::fopen(filename,"r");
    if (pin==NULL) std::cerr<<"Error: filename "<<filename<<" cannot be open!\n";
    else {
      step = NULL;
      bin_index = NULL;
      std::size_t rn;
      double a,b,g;
      rn = fread(&a, sizeof(double), 1 ,pin);
      if(rn<1) {
        std::cerr<<"Error: cannot read alpha!\n";
        abort();
      }
      rn = fread(&b, sizeof(double), 1 ,pin);
      if(rn<1) {
        std::cerr<<"Error: cannot read beta!\n";
        abort();
      }
      rn = fread(&g, sizeof(double), 1 ,pin);
      if(rn<1) {
        std::cerr<<"Error: cannot read gamma!\n";
        abort();
      }
      setabg(a,b,g);

      double error,dtm,dte;
      std::size_t itermax;
      int ext_method, ext_sequence;
      bool ext_iteration_const;
      
      rn = fread(&ext_method, sizeof(int), 1 ,pin);
      if(rn<1) {
        std::cerr<<"Error: cannot read exp_method!\n";
        abort();
      }
      rn = fread(&ext_sequence, sizeof(int), 1 ,pin);
      if(rn<1) {
        std::cerr<<"Error: cannot read exp_sequence!\n";
        abort();
      }
      rn = fread(&itermax, sizeof(std::size_t), 1 ,pin);
      if(rn<1) {
        std::cerr<<"Error: cannot read exp_itermax!\n";
        abort();
      }
      rn = fread(&error, sizeof(double), 1 ,pin);
      if(rn<1) {
        std::cerr<<"Error: cannot read exp_error!\n";
        abort();
      }
      rn = fread(&ext_iteration_const, sizeof(bool), 1 ,pin);
      if(rn<1) {
        std::cerr<<"Error: cannot read exp_fix_iter!\n";
        abort();
      }
      rn = fread(&dtm, sizeof(double), 1 ,pin);
      if(rn<1) {
        std::cerr<<"Error: cannot read dtmin!\n";
        abort();
      }
      rn = fread(&dte, sizeof(double), 1 ,pin);
      if(rn<1) {
        std::cerr<<"Error: cannot read dterr!\n";
        abort();
      }
      setEXP(error,dtm,dte,itermax,ext_method,ext_sequence,ext_iteration_const);

      int as;
      double asmin, asmax, aseps;
      std::size_t asitmin, asitmax;
      rn = fread(&as, sizeof(int), 1 ,pin);
      if(rn<1) {
        std::cerr<<"Error: cannot read auto_step!\n";
        abort();
      }
      rn = fread(&asmin, sizeof(double), 1 ,pin);
      if(rn<1) {
        std::cerr<<"Error: cannot read auto_step_fac_min!\n";
        abort();
      }
      rn = fread(&asmax, sizeof(double), 1 ,pin);
      if(rn<1) {
        std::cerr<<"Error: cannot read auto_step_fac_max!\n";
        abort();
      }
      rn = fread(&aseps, sizeof(double), 1 ,pin);
      if(rn<1) {
        std::cerr<<"Error: cannot read auto_step_eps!\n";
        abort();
      }
      rn = fread(&asitmin, sizeof(std::size_t), 1 ,pin);
      if(rn<1) {
        std::cerr<<"Error: cannot read auto_step_iter_min!\n";
        abort();
      }
      rn = fread(&asitmax, sizeof(std::size_t), 1 ,pin);
      if(rn<1) {
        std::cerr<<"Error: cannot read auto_step_iter_max!\n";
        abort();
      }
      setAutoStep(as,asmin,asmax,aseps,asitmin,asitmax);
    }
  }

  //! Print parameter table
  /*! Print parameter table
    @param[in] fout: ofstream for printing
   */
  void print(std::ostream & fout) {
    fout<<"Chain parameter table:\n"
        <<"Time step transformation parameters:\n"
        <<"  alpha = "<<alpha<<"  beta = "<<beta<<"  gamma = "<<gamma<<std::endl
        <<"Extrapolation parameters:\n"
        <<"  Interpolation method:    "<<((exp_method==1)?"Polynomial":"Rational")<<std::endl
        <<"  Sequence:                ";
    if (exp_sequence==1) fout<<"Romberg sequence {h, h/2, h/4, h/8 ...}\n";
    else if(exp_sequence==2) fout<<"Bulirsch & Stoer sequence {h, h/2, h/3, h/4, h/6, h/8 ...}\n";
    else if(exp_sequence==3) fout<<"4k sequence {h/2, h/6, h/10, h/14 ...}\n";
    else fout<<"Harmonic sequence {h, h/2, h/3, h/4 ...}\n";
    fout<<"  Maximum iteration times: "<<exp_itermax<<std::endl
        <<"  Phase/energy error criterion:    "<<exp_error<<std::endl
        <<"  Time sychronization error limit: "<<dterr<<std::endl
        <<"  Minimum physical time:           "<<dtmin<<std::endl
        <<"Auto-step parameters:\n"
        <<"  Auto-step method: ";
    if (auto_step==0) fout<<"None\n";
    else if(auto_step==1) fout<<"Use extrapolation error\n"
                              <<"  Minimum reduction factor: "<<auto_step_fac_min<<std::endl
                              <<"  Maximum increasing factor: "<<auto_step_fac_max<<std::endl;
    else if(auto_step==2) fout<<"Use min X/(g V) and V/(g A)\n"
                              <<"  Minimum reduction factor: "<<auto_step_fac_min<<std::endl
                              <<"  Maximum increasing factor: "<<auto_step_fac_max<<std::endl
                              <<"  Multiplied coefficient:   "<<auto_step_eps<<std::endl;
    else if(auto_step==3) fout<<"Use maximum sequence index to control next step\n"
                              <<"  Reduction factor: "<<auto_step_fac_min<<std::endl
                              <<"  Increasing factor: "<<auto_step_fac_max<<std::endl
                              <<"  Mimimum iteraction level: "<<auto_step_iter_min<<std::endl
                              <<"  Maximum iteraction level: "<<auto_step_iter_max<<std::endl;
    else if(auto_step==4) fout<<"Use minimum user-defined two-body timescale of each neigbor pairs\n"
                              <<"  Multiplied coefficient:   "<<auto_step_eps<<std::endl;
    fout<<"Iteration times fixed?: "<<exp_fix_iter<<std::endl;
  }
};

//! ARC class based on template class particle
/*!
  Major class for ARC integration of few bodies
  
  It depend on the template class particle. This particle class should contain public member functions for reading and writing mass, position and velocity (see sample in Particle::setPos(), Particle::setVel(), Particle::setMass(), Particle::getPos(), Particle::getVel(), Particle::getMass())
  The template depended type int_par is user-defined class for two-body interations used in ::ARC::chainpars.pair_AW, ::ARC::chainpars.pair_Ap


  The basic way to use ARC integration is shown as following:
  1. Construct a chain class with template class particle and a parameter controller of \ref ARC::chainpars. (The \ref ARC::chainpars should be configured first before doing integration. see its document for detail).
  2. Add existed particle 'A' (or a list of particles, or a chain type particle) into chain particle list (\ref chain.p) using chain.addP(). Notice the chain.addP() only registers the particle A's memory address into \ref chain.p without copying data. The chain integration will directly modify the position and velocity of particle A.
  3. Add perturbers into chain perturber list (\ref chain.pext) using chain.addPext() (also only register the particle address)
  4. Initialize Int_pars used for ::ARC::pair_AW or ::ARC::pair_Ap if necessary (see the case of ARC::Newtonian_AW())
  5. Initialize chain with chain.init(). Notice this function is necessary to be called before integration. Also be careful that after this initialization, the positions and velocites of particles registered in \ref chain.p will be shifted from their original frame to their center-of-mass frame. The particle type member variable \ref chain.cm stores the center-of-mass data of these particles (the mass of \ref chain.cm is the total mass of all member particles).
  6. Call integration functions (chain.Leapfrog_step_forward() or chain.extrapolation_integration()). The former use only Leapfrog method and the latter use extrapolation method to obtain high accuracy of integration.
  7. After call integration functions, the particles are integrated to new time. Because in ARC method, the time is also integrated and cannot be predicted before integration, thus the iteration need to be done to get correct physical time you want (see detailed in chain.extrapolation_integration() document).
  8. Notice that after chain.init() and integration, the particles are always in the center-of-mass frame. If you want to shift them back to the original frame, the chain.center_shift_inverse() should be used. But after this function is used, you should use chain.center_shift() before the next integration.

  Because of time now is an integrated variable, the time after integration cannot be predicted. 
  Thus if you want to stop the integration at a certain physical time, you need to use chain.extrapolation_integration() with dense output.
  To get better accuracy of physical time from intepolation, the 4k sequences (set in chainpars.setEXP()) is strongly suggested to be used. 
  If 4k sequences are used, the dense output method for GBS is used and the accuracy of time intepolation is close to the accuracy of integration.
  Please check the document of chain.extrapolation_integration() for detail.
 */
template <class particle, class int_par>
class chain{
  typedef double double3[3];
  double3 *X;  ///< relative position
  double3 *V;  ///< relative velocity
  double3 *acc; ///< acceleration
  double3 *pf;  ///< perturber force
  double3 *dWdr; ///< \partial Omega/ \partial rk
  std::size_t *list;   ///< chain index list

  //integration parameters=======================================//
  double t;    ///< time
  double Pt;    ///< Binding energy (time momentum)
  double w;    ///< integrated time transformation function

  //template parameters==========================================//
  double W;     ///< time transformation function
  double Ekin;  ///< kinetic energy
  double Pot;   ///< potential

  //number =======================================================//
  std::size_t num;      ///< total number of chain particles
  std::size_t nmax;     ///< maximum number 

  //monitor flags
  bool F_Pmod;     ///< indicate whether particle list is modified (true: modified)
  int  F_Porigin;  ///< indicate whether particle is shifted back to original frame (1: original frame: 0: center-of-mass frame; 2: only position is original frame)
  bool F_load;     ///< indicate whether load funcion is used

  const chainpars<int_par> *pars;   ///< chain parameter controller
  int_par* Int_pars;     ///< extra parameter array for pair_AW, pair_Ap and pair_T (used as last argument of ::ARC::pair_AW, ::ARC::pair_Ap and ::ARC::pair_T)

  chainlist<particle, int_par> p;    ///< particle list
  chainlist<particle, int_par> pext; ///< perturber list

public:

  particle cm;              ///< center mass particle

#ifdef TIME_PROFILE
  chainprofile profile;
#endif
  
  //! Constructor
  /*! Construct chain with allocated memory
      @param [in] n: maximum number of particles (will be used to allocate memory)
      @param [in] par: chain option controller class \ref ARC::chainpars
   */
  chain(const std::size_t n, const chainpars<int_par> &par):   num(0), pars(&par), Int_pars(NULL) {
    nmax=0;
    allocate(n);
  }

  //! Constructor
  /*! Construct chain without memory allocate, need to call allocate() later. 
     @param [in] par: chain option controller class \ref ARC::chainpars
   */
  chain(const chainpars<int_par> &par): num(0), nmax(0), F_Pmod(false), F_Porigin(1), F_load(false), pars(&par), Int_pars(NULL) {}

  //! Allocate memory
  /*! Allocate memory for maximum particle number n
     @param [in] n: maximum number of particles
   */
  void allocate(const std::size_t n) {
    if (nmax) {
      std::cerr<<"Error: chain memory allocation is already done\n";
      abort();
    }
    nmax = n;
    X=new double3[n-1];
    V=new double3[n-1];
    acc=new double3[n];
    pf=new double3[n];
    dWdr=new double3[n];
    list=new std::size_t[n];
    p.init(n);
    F_Pmod=false;
    F_Porigin=1;
    F_load=false;
  }

  //! Clear function
  /*! Clear allocated memory and set maximum number of particle to zero
   */
  void clear() {
    if (nmax>0) {
      delete[] X;
      delete[] V;
      delete[] acc;
      delete[] pf;
      delete[] dWdr;
      delete[] list;
      num = 0;
      nmax = 0;
    }
    p.clear();
    pext.clear();

    if (F_load) {
      if (Int_pars) delete Int_pars;
    }

    F_Pmod=false;
    F_Porigin=1;
    F_load=false;
#ifdef TIME_PROFILE
    profile.reset_tp();
#endif
  }

  //! destructor
  ~chain() {
    if (nmax>0) {
      delete[] X;
      delete[] V;
      delete[] list;
      delete[] acc;
      delete[] pf;
      delete[] dWdr;
    }
    if (F_load) {
      if (Int_pars) delete Int_pars;
    }

//    p.clear();
//    pext.clear();
  }

private:
//  //! Update number of particles
//  /*! Update the number of particle (#num) due to current particle list number
//     @param [in] n: current particle number in particle list #p
//  */
//  void update_num(const std::size_t n) {
//    if (n>nmax) {
//      std::cerr<<"Error: particle number "<<n<<" is larger than Chain number limit "<<num<<std::endl;
//      abort();
//    }
//    else{
//      num = n; ///< Update Chain number
//    }
//  }

  //! generate the chain list 
  /*! generate the chain list (#list) by N^2 searching the particle list (#p)
  */
  void generate_list() {
    bool *is_checked=new bool[num];
    for (std::size_t i=0; i<num; i++) is_checked[i] = false;
    std::size_t inext=0;
    for (std::size_t i=0; i<num; i++) {
      // initial rjk; mark checked particle
      is_checked[inext] = true;

      // initial chain_mem
      list[i]=inext;
      std::size_t inow=inext;
    
      // make chain
      double rmin;
      bool first=true;
      for (std::size_t j=1; j<num; j++) {
        if(is_checked[j]) continue;
        const double* rj = p[j].getPos();
        const double* ri = p[inow].getPos();
        double dx = rj[0] - ri[0];
        double dy = rj[1] - ri[1];
        double dz = rj[2] - ri[2];
        double dr2= dx*dx + dy*dy + dz*dz;
        if(first) {
          rmin = dr2;
          first=false;
          inext = j;
        }
        else if(dr2<rmin) {
          rmin = dr2;
          inext = j;
        }
      }
    }

    delete[] is_checked;
  }

  //! Calculate relative position and velocity
  /*! Get chain member relative position #X and velocity #V based on #list
  */
  void calc_XV() {
    for (std::size_t i=0;i<num-1;i++) {
      const double *ri1 = p[list[i+1]].getPos();
      const double *ri = p[list[i]].getPos();
      X[i][0] = ri1[0] - ri[0];
      X[i][1] = ri1[1] - ri[1];
      X[i][2] = ri1[2] - ri[2];

      const double *vi1 = p[list[i+1]].getVel();
      const double *vi = p[list[i]].getVel();
      V[i][0] = vi1[0] - vi[0];
      V[i][1] = vi1[1] - vi[1];
      V[i][2] = vi1[2] - vi[2];
    }
  }


  //! Calculate acceleration (potential) and transformation parameter
  /*! Get distance Matrix, acceleration, dm/dr, potential and transformation parameter
      based on particle masses in #p, using current #X, #V
      (notice the acceleration and dW/dr array index follow particle p to avoid additional shift when chain list change).

      @param [in] force: external force (acceleration) for each particle, (not perturber forces)
  */
  //      @param [in] resolve_flag: flag to determine whether to resolve sub-chain particles for force calculations. (defaulted false)
  void calc_rAPW (const double3 *force=NULL) {
#ifdef TIME_PROFILE
    profile.t_apw -= get_wtime();
#endif
    // safety check
    if (pars->pp_AW==NULL) {
      std::cerr<<"Error: acceleration and time transformation calculation function chainpars.pp_AW is not set!\n";
      abort();
    }
    // reset potential and transformation parameter
    double Pot_c  = 0.0;
    double W_c  = 0.0;
    // Loop all particles in list
#ifdef USE_OMP
#pragma omp parallel for reduction(+:Pot_c), reduction(+:W_c)
#endif
    for (std::size_t j=0;j<num;j++) {
      std::size_t lj = list[j];
      const particle *pj= &p[lj];

      for (std::size_t k=0; k<3; k++) {
        acc [lj][k]=0.0; // reset Acceleration
        dWdr[lj][k]=0.0; // reset dW/dr
      }

      for (std::size_t k=0;k<num;k++) {
        if(k==j) continue;
        std::size_t lk = list[k];
        const particle *pk= &p[lk];
        const double* xj = pj->getPos();
        double3 xjk;

        if(k==j+1) std::memcpy(xjk,X[j],3*sizeof(double));
        if(k==j-1) {
          xjk[0] = -X[k][0];
          xjk[1] = -X[k][1];
          xjk[2] = -X[k][2];
        }
        else if(k==j+2) {
          xjk[0] = X[j][0] + X[j+1][0];
          xjk[1] = X[j][1] + X[j+1][1];
          xjk[2] = X[j][2] + X[j+1][2];
        }
        else if(k==j-2) {
          xjk[0] = -X[k][0] - X[k+1][0];
          xjk[1] = -X[k][1] - X[k+1][1];
          xjk[2] = -X[k][2] - X[k+1][2];
        }
        else {
          const double* xk = pk->getPos();
          xjk[0] = xk[0] - xj[0];
          xjk[1] = xk[1] - xj[1];
          xjk[2] = xk[2] - xj[2];
        }

        double3 At,dWt;
        double Pt,Wt;
        const double mj=pj->getMass();
        const double mk=pk->getMass();

        // force calculation function from k to j
        pars->pp_AW(At, Pt, dWt, Wt, xjk, mj, mk, Int_pars);

//        // resolve sub-chain
//        if(resolve_flag && p.isChain(lk)) {
//          chain<particle, int_par>*ck = p.getSub(lk);
//          // center shift to current frame
//          ck->center_shift_inverse_X();
//          const std::size_t cn = ck->p.getN();
//          Pt = 0;
//          for (std::size_t i=0;i<3;i++) At[i]=0.0;
//          for (std::size_t i=0;i<cn;i++) {
//            double Ptemp;
//            double3 Atemp;
//            pars->pp_Ap(Atemp, Ptemp, xj, ck->p[i].getPos(), mj, ck->p[i].getMass());
// 
//            // Acceleration
//            At[0] += Atemp[0];
//            At[1] += Atemp[1];
//            At[2] += Atemp[2];
//            
//            // Potential
//            if (k>j) Pt += Ptemp;
//          }
//          // center shift back
//          ck->center_shift_X();
//        }

        // Acceleration
        acc[lj][0] += At[0];
        acc[lj][1] += At[1];
        acc[lj][2] += At[2];

        // dW/dr
        dWdr[lj][0] += dWt[0];
        dWdr[lj][1] += dWt[1];
        dWdr[lj][2] += dWt[2];

        if (k>j) {
          Pot_c += Pt; // Potential energy
          W_c += Wt;   // Transformation coefficient
        }
      }
      // add perturber force
      if (pext.getN()>0) {
        acc[lj][0] += pf[lj][0];
        acc[lj][1] += pf[lj][1];
        acc[lj][2] += pf[lj][2];
      }        
      // add external acceleration
      if (force!=NULL) {
        acc[lj][0] += force[lj][0];
        acc[lj][1] += force[lj][1];
        acc[lj][2] += force[lj][2];
      }
    }

    // Write potential and w
    Pot = Pot_c;
    W = W_c;
#ifdef TIME_PROFILE
    profile.t_apw += get_wtime();
#endif
  }

  //! Calculate kinetic energy
  void calc_Ekin(){
    Ekin = 0.0;
    for (std::size_t i=0; i<num; i++) {
      const double *vi=p[i].getVel();
      Ekin += 0.5 * p[i].getMass() * (vi[0]*vi[0]+vi[1]*vi[1]+vi[2]*vi[2]);
    }
  }

  //! Step forward of X
  /*! One step integration of #X 
     @param [in] dt: physical time step dt for #X
   */
  void step_forward_X(const double dt) {
    // step forward relative X
    for (std::size_t i=0;i<num-1;i++) {
      X[i][0] += dt * V[i][0];
      X[i][1] += dt * V[i][1];
      X[i][2] += dt * V[i][2];
    }
  }

  //! Step forward of V
  /*! one step integration of V
     @param [in] dt: physical time step dt for #V
   */
  void step_forward_V(const double dt) {
    // step forward V
    for (std::size_t i=0;i<num-1;i++) {
      std::size_t k = list[i];
      std::size_t k1 = list[i+1];
      V[i][0] += dt * (acc[k1][0]-acc[k][0]);
      V[i][1] += dt * (acc[k1][1]-acc[k][1]);
      V[i][2] += dt * (acc[k1][2]-acc[k][2]);
    }
  }

  //! Step forward of Pt and w
  /*! One step integration of Pt and w.
      - \f$Pt += dt * \sum ( - m_k * <v_k> \dot f_k)\f$ 
      - \f$w += dt * \sum ( dm/dr_k \dot v_k)\f$ 
     @param [in] dt: time step for V
     @param [in] ave_v: averaged velocity
     @param [in] force: external force
     @param [in] p: particle list (only use mass)
     @param [in] fpf: if perturber force is not zero, it is true
   */
  void step_forward_Ptw(const double dt, const double3* ave_v, const double3* force, const bool fpf) {
    double dPt = 0.0;
    double dw = 0.0;
    if (force!=NULL||fpf||pars->beta>0) {
      for (std::size_t i=0;i<num;i++) {
        if (force!=NULL) {
          dPt -= p[i].getMass() * ( ave_v[i][0] * (pf[i][0] + force[i][0]) 
                            + ave_v[i][1] * (pf[i][1] + force[i][1]) 
                            + ave_v[i][2] * (pf[i][2] + force[i][2]));
        }
        else if (fpf){
          dPt -= p[i].getMass() * ( ave_v[i][0] * pf[i][0] 
                            + ave_v[i][1] * pf[i][1] 
                            + ave_v[i][2] * pf[i][2]);
        }
        
        if (pars->beta>0) {
          dw += ( ave_v[i][0] * dWdr[i][0]
                + ave_v[i][1] * dWdr[i][1]
                + ave_v[i][2] * dWdr[i][2]);
        }
      }
    }
    Pt += dt * dPt;
    w += dt * dw;
  }

  //! resolve X and V
  /*! resolve relative #X, #V to physical x, v and calculated the averaged velocity of old and new values.
      Notice the center-of-mass particle mass in Chain.cm is used.
      The total mass of particles should be consistent with cm.getMass(). Otherwise update Chain.cm first.
      @param [out] ave_v: averaged velocity array (return values)
   */
  void resolve_XV(double3* ave_v=NULL) {
    // backup old v
    if (ave_v!=NULL) {
      for (std::size_t i=0;i<num;i++) {
        const double *vi = p[i].getVel();
        ave_v[i][0] = vi[0];
        ave_v[i][1] = vi[1];
        ave_v[i][2] = vi[2];
      }
    }
    // resolve current X,V
    // first member
    double3 vc;
    double3 xc;
    const std::size_t lk1=list[0];
    const double  mk1 = p[lk1].getMass();
    const double *rk1 = p[lk1].getPos();
    const double *vk1 = p[lk1].getVel();
    xc[0] = mk1 * rk1[0];
    xc[1] = mk1 * rk1[1];
    xc[2] = mk1 * rk1[2];
    vc[0] = mk1 * vk1[0];
    vc[1] = mk1 * vk1[1];
    vc[2] = mk1 * vk1[2];
    // others
    for (std::size_t i=0;i<num-1;i++) {
      const std::size_t lk = list[i];
      const std::size_t lkn = list[i+1];
      const double *rk = p[lk].getPos();
      p[lkn].setPos(rk[0] + X[i][0],
                    rk[1] + X[i][1],
                    rk[2] + X[i][2]);

      const double *vk = p[lk].getVel();
      p[lkn].setVel(vk[0] + V[i][0],
                    vk[1] + V[i][1],
                    vk[2] + V[i][2]);

      //! center-of-mass position and velocity
      const double mkn = p[lkn].getMass();
      const double *rkn = p[lkn].getPos();
      const double *vkn = p[lkn].getVel();
      xc[0] += mkn * rkn[0];
      xc[1] += mkn * rkn[1];
      xc[2] += mkn * rkn[2];
      vc[0] += mkn * vkn[0];
      vc[1] += mkn * vkn[1];
      vc[2] += mkn * vkn[2];
    }

    // calcualte center-of-mass position and velocity shift
    const double mc = cm.getMass();
    xc[0] /= mc;
    xc[1] /= mc;
    xc[2] /= mc;
    vc[0] /= mc;
    vc[1] /= mc;
    vc[2] /= mc;
    
    for (std::size_t i=0;i<num;i++) {
      // center-of-mass correction
      const double *ri = p[i].getPos();
      p[i].setPos(ri[0] - xc[0],
                  ri[1] - xc[1],
                  ri[2] - xc[2]);
      const double *vi = p[i].getVel();
      p[i].setVel(vi[0] - vc[0],
                  vi[1] - vc[1],
                  vi[2] - vc[2]);

      // calculate averaged velocities
      if (ave_v!=NULL) {
        ave_v[i][0] = 0.5 * (ave_v[i][0] + vi[0]);
        ave_v[i][1] = 0.5 * (ave_v[i][1] + vi[1]);
        ave_v[i][2] = 0.5 * (ave_v[i][2] + vi[2]);
      }
    }
  }

  //! Resolve X
  /*! Resolve relative #X to physical x (center-of-mass frame)
      Notice the center-of-mass particle mass in #Chain.cm is used.
      The total mass of particles should be consistent with cm.getMass(). Otherwise update Chain.cm first.
   */
  void resolve_X() {
    // resolve current X
    double3 xc;
    // first member
    const std::size_t lk1=list[0];
    const double  mk1 = p[lk1].getMass();
    const double *rk1 = p[lk1].getPos();
    xc[0] = mk1 * rk1[0];
    xc[1] = mk1 * rk1[1];
    xc[2] = mk1 * rk1[2];
    // others
    for (std::size_t i=0;i<num-1;i++) {
      const std::size_t lk = list[i];
      const std::size_t lkn = list[i+1];
      const double *rk = p[lk].getPos();
      p[lkn].setPos(rk[0] + X[i][0],
                    rk[1] + X[i][1],
                    rk[2] + X[i][2]);

      // center-of-mass position and velocity
      const double mkn = p[lkn].getMass();
      const double *rkn = p[lkn].getPos();
      xc[0] += mkn * rkn[0];
      xc[1] += mkn * rkn[1];
      xc[2] += mkn * rkn[2];
    }

    // calcualte center-of-mass position and velocity shift
    const double mc = cm.getMass();
    xc[0] /= mc;
    xc[1] /= mc;
    xc[2] /= mc;
    
    for (std::size_t i=0;i<num;i++) {
      // center-of-mass correction
      const double *ri = p[i].getPos();
      p[i].setPos(ri[0] - xc[0],
                  ri[1] - xc[1],
                  ri[2] - xc[2]);
    }
  }

  //! Resolve V
  /*! Resolve relative velocity #V to physical  v (center-of-mass frame)
    Notice the center-of-mass particle mass in Chain.cm is used.
    The total mass of particles should be consistent with cm.getMass(). Otherwise update Chain.cm first.
   */
  void resolve_V() {
    // resolve current V
    // first member
    double3 vc;
    const std::size_t lk1=list[0];
    const double  mk1 = p[lk1].getMass();
    const double *vk1 = p[lk1].getVel();
    vc[0] = mk1 * vk1[0];
    vc[1] = mk1 * vk1[1];
    vc[2] = mk1 * vk1[2];
    // others
    for (std::size_t i=0;i<num-1;i++) {
      const std::size_t lk = list[i];
      const std::size_t lkn = list[i+1];
      const double *vk = p[lk].getVel();
      p[lkn].setVel(vk[0] + V[i][0],
                    vk[1] + V[i][1],
                    vk[2] + V[i][2]);

      // center-of-mass position and velocity
      const double mkn = p[lkn].getMass();
      const double *vkn = p[lkn].getVel();
      vc[0] += mkn * vkn[0];
      vc[1] += mkn * vkn[1];
      vc[2] += mkn * vkn[2];
    }

    // calcualte center-of-mass position and velocity shift
    const double mc = cm.getMass();
    vc[0] /= mc;
    vc[1] /= mc;
    vc[2] /= mc;
    
    for (std::size_t i=0;i<num;i++) {
      // center-of-mass correction
      const double *vi = p[i].getVel();
      p[i].setVel(vi[0] - vc[0],
                  vi[1] - vc[1],
                  vi[2] - vc[2]);
    }
  }

  //! Perturber force calculation
  /*! Get perturber force based on porturber list #pext
    @param [in] resolve_flag: whether resolve perturber member if it is a chain 
    \return flag: true: pertubers exist. false: no perturbers
  */
  bool pert_force(const bool resolve_flag=false) {
#ifdef TIME_PROFILE
    profile.t_pext -= get_wtime();
#endif
    // safety check
    if (pars->pp_Ap==NULL) {
      std::cerr<<"Error: acceleration calculation function chainpars.pp_Ap is not set!\n";
      abort();
    }
    const int np = pext.getN();
    if (np>0) {
      for (std::size_t i=0;i<num;i++) {
        for (std::size_t j=0;j<3;j++) pf[i][j] = 0.0;
        for (std::size_t j=0;j<np;j++) {
          double3 At={0};
          double Pt;
          const particle* pi=&p[i];
          const double* xi=pi->getPos();
          const double  mi=pi->getMass();
          
          // check sub-chain system
          if (resolve_flag && pext.isChain(j)) {
            chain<particle, int_par>*cj = pext.getSub(j);
            // get center-of-mass position for shifting;
            const double* xc=cj->cm.getPos();
            const std::size_t cn = cj->pext.getN();
            for (std::size_t k=0;k<cn;k++) {

              double3 xk;
              std::memcpy(xk,cj->p[k].getPos(),3*sizeof(double));
              // shift position to current frame; to keep thread safety, original data are not modified
              if (cj->isPorigin()==0) {
                xk[0] += xc[0];
                xk[1] += xc[1];
                xk[2] += xc[2];
              }
              double3 Atemp;
              pars->pp_Ap(Atemp, Pt, xi, xk, mi, cj->p[k].getMass(), Int_pars);

              // Acceleration
              At[0] += Atemp[0];
              At[1] += Atemp[1];
              At[2] += Atemp[2];
            }
          }
          else {
            // perturber force
            pars->pp_Ap(At, Pt, xi, pext[j].getPos(), mi, pext[j].getMass(), Int_pars);
          }
          
          pf[i][0] += At[0];
          pf[i][1] += At[1];
          pf[i][2] += At[2];
        }
      }
#ifdef TIME_PROFILE
      profile.t_pext += get_wtime();
#endif
      return true;
    }
    else {
      std::memset(pf,0,3*num*sizeof(double));
#ifdef TIME_PROFILE
      profile.t_pext += get_wtime();
#endif
      return false;
    }
  }
     
  //! Update #list order based on the relative distances
  /*! Update chain #list order based on current relative distance #X
    \return flag: if link is modified, return true
   */
  bool update_link(){
#ifdef TIME_PROFILE
    profile.t_uplink -= get_wtime();
#endif
    bool modified=false; // indicator
#ifdef DEBUG        
    std::cerr<<"current:";
    for (std::size_t i=0;i<num;i++) std::cerr<<std::setw(4)<<list[i];
    std::cerr<<"\n";
#endif
    
    // create reverse index of link
    std::size_t* rlink = new std::size_t[num];
    std::size_t* roldlink = new std::size_t[num];
    for (std::size_t i=0;i<num;i++) rlink[list[i]] = i;
    std::memcpy(roldlink,rlink,num*sizeof(std::size_t));

    // backup previous link
    std::size_t* listbk = new std::size_t[num];
    std::memcpy(listbk,list,num*sizeof(std::size_t));

    // backup current X
    double3* Xbk = new double3[num-1];
    std::memcpy(Xbk[0],X[0],(num-1)*3*sizeof(double));
    
    // backup current V
    double3* Vbk = new double3[num-1];
    std::memcpy(Vbk[0],V[0],(num-1)*3*sizeof(double));

    // create mask to avoid dup. check;
    bool* mask = new bool[num];
    for (std::size_t i=0;i<num;i++) mask[i] = false;

    const double NUMERIC_DOUBLE_MAX = std::numeric_limits<double>::max();
    for (std::size_t k=0;k<num-1;k++) {
      std::size_t lk  = list[k];
       mask[lk] = true;
      std::size_t lkn = list[k+1];
      // possible new index
      std::size_t lku = lkn;
      // calculate distance
      double rmin = NUMERIC_DOUBLE_MAX;
      for (std::size_t j=0;j<num;j++) {
        if (mask[j]||j==k) continue;
        const double* xk = p[lk].getPos();
        const double* xj = p[j].getPos();
        double xjk = xk[0] - xj[0];
        double yjk = xk[1] - xj[1];
        double zjk = xk[2] - xj[2];
        double rjk2 = xjk*xjk + yjk*yjk + zjk*zjk;
        if (rjk2<rmin) {
          lku=j;
          rmin = rjk2;
        }
      }
      if (lku!=lkn) {
#ifdef DEBUG        
        std::cerr<<"Switch: "<<k<<" new: "<<lku<<" old: "<<lkn<<std::endl;
        for (std::size_t i=0;i<num;i++) std::cerr<<std::setw(4)<<list[i];
        std::cerr<<"\n Rlink";
        for (std::size_t i=0;i<num;i++) std::cerr<<std::setw(4)<<rlink[i];
        std::cerr<<"\n";
#endif
        modified=true;
        // shift two index in the list
        list[rlink[lku]] = lkn;
        rlink[lkn] = rlink[lku];
        list[k+1] = lku;
        rlink[lku] = k+1;
        mask[lku] = true;
#ifdef DEBUG        
        for (std::size_t i=0;i<num;i++) std::cerr<<std::setw(4)<<list[i];
        std::cerr<<"\n Rlink";
        for (std::size_t i=0;i<num;i++) std::cerr<<std::setw(4)<<rlink[i];
        std::cerr<<"\n";
#endif
      }

      if (lk!=listbk[k]||lku!=listbk[k+1]) {
        // update X and V
        // left boundary
        std::size_t rlk = roldlink[lk];
        if (rlk<k) {
          for (std::size_t j=rlk;j<k;j++) {
#ifdef DEBUG
            std::cerr<<"Add V["<<j<<"] to V["<<k<<"]\n";
#endif
            X[k][0] += Xbk[j][0];
            X[k][1] += Xbk[j][1];
            X[k][2] += Xbk[j][2];

            V[k][0] += Vbk[j][0];
            V[k][1] += Vbk[j][1];
            V[k][2] += Vbk[j][2];
          }
        }
        else if (rlk>k) {
          for (std::size_t j=k;j<rlk;j++) {
#ifdef DEBUG
            std::cerr<<"Minus V["<<j<<"] from V["<<k<<"]\n";
#endif
            X[k][0] -= Xbk[j][0];
            X[k][1] -= Xbk[j][1];
            X[k][2] -= Xbk[j][2];

            V[k][0] -= Vbk[j][0];
            V[k][1] -= Vbk[j][1];
            V[k][2] -= Vbk[j][2];
          }
        }
        // right boundary
        rlk = roldlink[lku];
        if (rlk<k+1) {
          for (std::size_t j=rlk;j<k+1;j++) {
#ifdef DEBUG
            std::cerr<<"Minus V["<<j<<"] from V["<<k<<"]\n";
#endif
            X[k][0] -= Xbk[j][0];
            X[k][1] -= Xbk[j][1];
            X[k][2] -= Xbk[j][2];

            V[k][0] -= Vbk[j][0];
            V[k][1] -= Vbk[j][1];
            V[k][2] -= Vbk[j][2];
          }
        }
        else if (rlk>k+1) {
          for (std::size_t j=k+1;j<rlk;j++) {
#ifdef DEBUG
            std::cerr<<"Add V["<<j<<"] to V["<<k<<"]\n";
#endif
            X[k][0] += Xbk[j][0];
            X[k][1] += Xbk[j][1];
            X[k][2] += Xbk[j][2];

            V[k][0] += Vbk[j][0];
            V[k][1] += Vbk[j][1];
            V[k][2] += Vbk[j][2];
          }
        }
      }
    }

    // clear template array
    delete[] mask;
    delete[] rlink;
    delete[] roldlink;
    delete[] listbk;
    delete[] Vbk;
    delete[] Xbk;

#ifdef DEBUG    
    if(modified) print();
#endif
    
#ifdef TIME_PROFILE
    profile.t_uplink += get_wtime();
#endif
    return modified;
  }

  //! Center of mass shift 
  /*! Shift positions and velocities of N (#num) particles (#p) based on their center-of-mass, write center-of-mass particle to chain
  */
  void center_shift_init() {
    // center mass
    double3 cmr={};
    double3 cmv={};
    double  cmm=0;
    for (std::size_t i=0;i<num;i++) {
      const double *ri = p[i].getPos();
      const double *vi = p[i].getVel();
      const double mi = p[i].getMass();
      cmr[0] += ri[0] * mi;
      cmr[1] += ri[1] * mi;
      cmr[2] += ri[2] * mi;

      cmv[0] += vi[0] * mi;
      cmv[1] += vi[1] * mi;
      cmv[2] += vi[2] * mi;

      cmm += mi;
    }
    cmr[0] /= cmm; 
    cmr[1] /= cmm; 
    cmr[2] /= cmm; 

    cmv[0] /= cmm; 
    cmv[1] /= cmm; 
    cmv[2] /= cmm;

    cm.setMass(cmm);
    cm.setPos(cmr[0],cmr[1],cmr[2]);
    cm.setVel(cmv[0],cmv[1],cmv[2]);

    // shifting
    for (std::size_t i=0;i<num;i++) {
      const double *ri = p[i].getPos();
      const double *vi = p[i].getVel();
      p[i].setPos(ri[0] - cmr[0],
                  ri[1] - cmr[1],
                  ri[2] - cmr[2]);
      p[i].setVel(vi[0] - cmv[0],
                  vi[1] - cmv[1],
                  vi[2] - cmv[2]);
    }
  }

  //! Inversed center of mass shift for #X
  /*! Shift back positions of N (num) particles (p) based on chain center-of-mass
  */
  void center_shift_inverse_X() {
    if (F_Porigin>0) {
      std::cerr<<"Warning: particles are already in original frame!\n";
    }
    else {
      const double *rc = cm.getPos();
      for (std::size_t i=0;i<num;i++) {
        const double *ri = p[i].getPos();
        p[i].setPos(ri[0] + rc[0],
                    ri[1] + rc[1],
                    ri[2] + rc[2]);
      }
      F_Porigin=2;
    }
  }

  //! Center of mass shift for #X
  /*! Shift positions and velocities of N (num) particles (p) based on chain center-of-mass
  */
  void center_shift_X() {
    if (F_Porigin==0) {
      std::cerr<<"Warning: particles are already in the center-of-mass frame!\n";
    }
    else if (F_Porigin==1) {
      std::cerr<<"Error: particles are in original frame, please use center_shift instead!\n";
      abort();
    }
    else {
      const double *rc = cm.getPos();
      for (std::size_t i=0;i<num;i++) {
        const double *ri = p[i].getPos();
        p[i].setPos(ri[0] - rc[0],
                    ri[1] - rc[1],
                    ri[2] - rc[2]);
      }
      F_Porigin=0;
    }
  }


  //! Middle difference calculator
  /*! Central difference cumulative calculator for time #t
     If \a i is 0, the initial values will be reset to zero, and for the same n, this function need to be called n+1 times to complete the n'th order difference
    @param [in] dpoly: two dimensional array to storage the difference of dataset for different order. Array size is [\a nmax][dsize]. dsize = 1 due to one point
                       If dpoly is NULL, no calculation will be done.
    @param [in] nmax: maximum difference level (count from 1)
    @param [in] i: current position of dataset corresponding to the binomial coefficients (count from 0 as starting point)
    @param [in] ndiv: substep number
  */
  void mid_diff_calc(double **dpoly, const std::size_t nmax, const std::size_t i, const int ndiv) {
    // safety check
    if (dpoly!=NULL) {
      // difference level should not exceed the point number
      if (2*nmax>ndiv) {
        std::cerr<<"Error: maximum difference order "<<nmax<<" *2 > total step number "<<ndiv<<"!\n";
        abort();
      }
      // for middle difference, odd substep number is not allown
      if(ndiv%2) {
        std::cerr<<"Error: middle difference cannot allow odd substep number ("<<ndiv<<")!\n";
        abort();
      }
      if (i>=ndiv/2-nmax&&i<=ndiv/2+nmax) {
//        if (!tflag) {
//          const std::size_t dsize=6*num-3;
//          //initial dpoly to zero
//          for (std::size_t j=0; j<nmax; j++) if(i==0) std::memset(dpoly[j],0,dsize*sizeof(double));
//        }
//        else
        if(i==0) for (std::size_t j=0; j<nmax; j++) dpoly[j][0] = 0.0;

        // dt/ds
        double dts=calc_dt_X(1.0);
        
        int** binI = pars->bin_index;
        // formula delta^n f(x) = sum_ik=0,n (-1)^(n-ik) (n n-ik) f(x+2*ik*h)
        for (int j=0; j<nmax; j++) {
          if ((i%2&&j%2)||(i%2==0&&j%2==0)) {
            // n=j+1 indicate the difference degree (first difference is j=0)
            const int n=j+1;
            // for even and odd difference level n, different i points are used
            // ndiv/2: middle point; n/2: left edge point shift, (should be double since only odd points are used);
            int ik= i-(ndiv/2-n);
            // ik should be limited for diferent level of difference
            if (ik>=0&&ik<=2*n) {
              int nk= n-ik/2;  // n-ik
              double coff = ((nk%2)?-1:1)*binI[n][nk];
//#ifdef DEBUG
//              std::cerr<<"Poly_coff: n= "<<n<<"; i="<<i<<"; ik="<<ik<<"; coff="<<coff<<std::endl;
//#endif
              dpoly[j][0] += coff * dts;
//              if (!tflag) {
//                dpoly[j][1] += coff * Pt;
//                dpoly[j][2] += coff * w;
//                for (std::size_t k=0; k<num-1; k++) {
//                  for (std::size_t kk=0; kk<3; kk++) {
//                    dpoly[j][3*(1+k)+kk] += coff * X[k][kk];
//                    dpoly[j][3*(num+k)+kk] += coff * V[k][kk];
//                  }
//                }
//              }
            }
          }
        }
      }
    }
  }


  //! Edge difference calculator
  /*! Left forward and right backward difference cumulative calculator for time #t
      If i is 0, the array will be reset to zero, and for the same n, this function need to be called n+1 times (i=0...n+1) to complete the n'th order difference.
      @param [out] dpoly: two dimensional array to storage the results. Array size is [\a nmax][dsize].  dsize =2 because of two points.
                          first [] indicate different level of differences, second [] indicate the left and right difference of time #t
                          If dpoly is NULL, no calculation will be done.    
      @param [in] nmax: maximum difference level (count from 1)
      @param [in] i: current position of data corresponding to the binomial coefficients (count from 0 as starting point)              
      @param [in] ndiv: substep number
t)
   */
  void edge_diff_calc(double **dpoly, const int nmax, const std::size_t i, const int ndiv) {
    if (dpoly!=NULL) {
      // safety check
//      if (!tflag) {
//        const std::size_t dsize=12*num-6;
//        // i indicate the position, i=0 means x0
//        if(i==0) for (int j=0; j<nmax; j++) std::memset(dpoly[j],0,dsize*sizeof(double)); // reset data array at i=0
//      }
//      else
      if(i==0) for (int j=0; j<nmax; j++) {
          dpoly[j][0]=0.0;
          dpoly[j][1]=0.0;
        }
      
      if (nmax>ndiv) {
        std::cerr<<"Error: maximum difference order "<<nmax<<" > total step number "<<ndiv<<"!\n";
        abort();
      }

      if (i<=nmax||i>=ndiv-nmax) {
        // dt/ds
        double dts=calc_dt_X(1.0);

        int** binI = pars->bin_index;
        for (int j=0; j<nmax; j++) {
          // j+1 indicate the difference degree, count from 1 (first difference)
          const int n = j+1;
              
          // left edge, forward difference: (-1)^(n-i) (n n-i) f(x0+i*h) (i from 0 to n)
          int ileft = (int)(n-i);   // n-i
          if (ileft>=0) {
            double coff = ((ileft%2)?-1:1)*binI[n][ileft];
            dpoly[j][0] += coff * dts;
//            if (!tflag) {
//              dpoly[j][1] += coff * Pt;
//              dpoly[j][2] += coff * w;
//              for (std::size_t k=0; k<num-1; k++) {
//                for (std::size_t kk=0; kk<3; kk++) {
//                  dpoly[j][3*(1+k)+kk] += coff * X[k][kk];
//                  dpoly[j][3*(num+k)+kk] += coff * V[k][kk];
//                }
//              }
//            }
#ifdef DEBUG
            std::cerr<<"Poly left: n="<<n<<" ik="<<ileft<<" i="<<i<<" coff="<<coff<<" t="<<t<<std::endl;
#endif
          }

          // right edge, backward difference: (-1)^(n-i) (n n-i) f(xn-(n-i)*h) (i from 0 to n)
          int ishift = ndiv-n;
          int iright= ileft+ishift;     // n-i
          if (i>=ishift) {
            double coff = ((iright%2)?-1:1)*binI[n][iright];
//            if (!tflag) {
//              const std::size_t ir = 6*num-3;
//              dpoly[j][0+ir] += coff * t;
//              dpoly[j][1+ir] += coff * Pt;
//              dpoly[j][2+ir] += coff * w;
//              for (std::size_t k=0; k<num-1; k++) {
//                for (std::size_t kk=0; kk<3; kk++) {
//                  dpoly[j][3*(1+k)+kk+ir] += coff * X[k][kk];
//                  dpoly[j][3*(num+k)+kk+ir] += coff * V[k][kk];
//                }
//              }
//            } else
            dpoly[j][1] += coff * dts;
#ifdef DEBUG
            std::cerr<<"                 Poly right: n="<<n<<" ik="<<iright<<" i="<<i<<" coff="<<coff<<std::endl;
#endif
          }
        }
      }
    }
  }

  //! diff_dev_calc
  /*£¡ Calcualte derivates from differences: \d^(n)f/ h^n
    @param [in,out] dpoly: two dimensional storing differences, will be updated to derivates. Array size is [nmax][dsize]. 
                          first [] indicate different level of differences, second [] indicate the difference of data. 
                          If dpoly is NULL, no calculation will be done.
    @param [in] h: step size
    @param [in] nmax: maximum difference order
    @param [in] dsize: number of data in dataset
   */
  void diff_dev_calc(double **dpoly, const double h, const int nmax, const int dsize) {
    double hn = h;
    // loop difference order from 1 to nmax
    for (std::size_t i=0; i<nmax; i++) {
#ifdef DEBUG
      std::cerr<<"Diff order "<<i<<"; h="<<hn<<"; nmax="<<nmax<<std::endl;
#endif
      for (std::size_t j=0; j<dsize; j++) dpoly[i][j] /= hn;
      hn *= h;
    }
  }

public:
  //! Find the most strongly interacted pair
  /*! Find the most strongly interacted pair.
    First the dacc[i] = acc[i+1]-acc[i] (i=0, num-1) is calculated, then select the maximum dacc[i], the corresponding two particles are selected as the strong pair. Their indice are written to pindex.
    @param [in] pindex: two element array used for storing pair index.
    \return the acceleraction difference of this pair
   */
  double find_strong_pair(int pindex[2]) {
    std::size_t k=0;
    double accm = 0;
    for (std::size_t i=0; i<num-1; i++) {
      const double accx = acc[i+1][0]-acc[i][0];
      const double accy = acc[i+1][1]-acc[i][1];
      const double accz = acc[i+1][2]-acc[i][2];
      double acci = accx*accx + accy*accy* + accz*accz;
      if (accm<acci) {
        accm = acci;
        k = i;
      }
    }
    pindex[0] = list[k];
    pindex[1] = list[k+1];

    return std::sqrt(accm);
  }

  //! Find the closest pair
  /*! Find the closest pair using relative positions
    @param [in] pindex: two element array used for storing pair index.
    \return the relative distance of this pair
   */
  double find_closest_pair(int pindex[2]) {
    std::size_t k=0;
    double Xm = std::numeric_limits<double>::max();
    for (std::size_t i=0; i<num-1; i++) {
      double Xi2 = X[i][0]*X[i][0] + X[i][1]*X[i][1] + X[i][2]*X[i][2];
      if (Xm > Xi2) {
        Xm = Xi2;
        k = i;
      }
    }
    pindex[0] = list[k];
    pindex[1] = list[k+1];

    return std::sqrt(Xm);
  }
  
  
  //! Calculate physical time step for X
  /*! Calculate physical time step dt for #X based on ds
     @param [in] ds:  integration step size (not physical time step) 
     \return     dt:  physical integration time step for #X
  */
  double calc_dt_X(const double ds) {
    // determine the physical time step
    double dt = ds / (pars->alpha * (Ekin + Pt) + pars->beta * w + pars->gamma);
    if (std::abs(dt) < pars->dtmin) {
      std::cerr<<"Error!: physical time step too small: "<<dt<<std::endl;
      abort();
    }
    return dt;
  }
  
  //! Calculate physical time step for V
  /*! Calculate physical time step dt for #V based on ds
     @param [in] ds:  step size s (not physical time step) 
     \return     dt:  physical integration time step for #V
  */
  double calc_dt_V(const double ds) {
    // determine velocity integration time step
    return ds / (pars->gamma - pars->alpha * Pot + pars->beta * W);
  }
  

  //! Calculate next step approximation based on min(X/(gV),V/(gA))
  /*!
    \return: approximation of step size ds
  */
  double calc_next_step_XVA() {
    double dsXV=std::numeric_limits<double>::max();
    double dsVA=std::numeric_limits<double>::max();
    for (std::size_t i=0; i<num-1; i++) {
      double r = X[i][0]*X[i][0]+X[i][1]*X[i][1]+X[i][2]*X[i][2];
      double v = V[i][0]*V[i][0]+V[i][1]*V[i][1]+V[i][2]*V[i][2];
      double a = acc[i][0]*acc[i][0]+acc[i][1]*acc[i][1]+acc[i][2]*acc[i][2];
      dsXV = std::min(r/v,dsXV);
      dsVA = std::min(v/a,dsVA);
    }
    return pars->auto_step_eps*std::min(std::sqrt(dsXV)/calc_dt_X(1.0),std::sqrt(dsVA)/calc_dt_V(1.0));
  }

  //! Calculate next step approximation based on custom defined timescale for two-body system
  /*!
    The user-defined two-body timescale function with ::ARC::pair_T type will be performed on all neighbor pairs in chain.
    Then the minimum timescale \f$ T_m\f$ is selected to calculate next step size: \f$ ds = \eps T_m |Pt|\f$, where \f$\eps\f$ is #auto_step_eps set in chainpars.setAutoStep().
    @param [in] pt: two-body timescale function with ::ARC::pair_T type
    \return approximation of step size ds
   */
  double calc_next_step_custom() {
    //safety check
    if (pars->pp_T==NULL) {
      std::cerr<<"Error: two-body timescale calculator function chainpars.pp_T is not set\n";
      abort();
    }

    double perim = std::numeric_limits<double>::max();
    for (std::size_t i=0; i<num-1; i++) {
      const double peri=pars->pp_T(p[list[i]].getMass(),p[list[i+1]].getMass(),X[i],V[i],Int_pars);
      if (perim>peri) perim = peri;
    }
    return pars->auto_step_eps*perim*std::abs(Pt);
  }
  
  //! Add particle
  /*! Add one particle (address pointer) into particle list #p (see ARC::chainlist.add())
     @param [in] a: new particle
   */
  void addP(particle &a) {
    if (F_Porigin!=1) std::cerr<<"Warning!: particle list are (partically) in the center-of-mass frame, dangerous to add new particles!\n";
    p.add(a);
    F_Pmod=true;
  }
  
  //! Add chain as a particle in #p
  /*! Add one chain (address pointer) into particle list #p (see ARC::chainlist.add())
    @param [in] a: new chain particle
   */
  void addP(chain<particle,int_par> &a) {
    if (F_Porigin!=1) std::cerr<<"Warning!: particle list are (partically) in the center-of-mass frame, dangerous to add new particles!\n";
    p.add(a);
    F_Pmod=true;
  }
  
  //! Add a list of particle
  /*! Add a list of particles (see ARC::chainlist.add())
    @param [in] n: number of particles need to be added
    @param [in] a: array of new particles
   */
  void addP(const std::size_t n, particle a[]) {
    if (F_Porigin!=1) std::cerr<<"Warning!: particle list are (partically) in the center-of-mass frame, dangerous to add new particles!\n";
    p.add(n,a);
    F_Pmod=true;
  }

  //! remove one particle
  /*! remove one particle from #p (see ARC::chainlist.remove())
     @param [in] i: particle index in #p needs to be removed
     @param [in] option: position update option: 
                 - true: shift last particle to current position (defaulted);
                 - false: shift all right particle to left by one
  */
  void removeP(const std::size_t i, bool option=true) { p.remove(i,option); F_Pmod=true; }

  //! Data Dumping to file function
  /*! Dump chain data into file (binary format) using fwrite. 
      The first variable in the dumping data is particle number #num
      The second is the particle chain list index array #list
      The third is one dimensional array data storing #t, #B, #w, #X, #V generated by using backup().
      The forth is center of mass particle data
      The fifth is one dimensional array storing the masses of particles
      @param [in] filename: file for storing the data
   */
  void dump(const char* filename) {
    std::FILE* pout = std::fopen(filename,"w");
    if (pout==NULL) std::cerr<<"Error: filename "<<filename<<" cannot be open!\n";
    else {
      // number of particles
      fwrite(&num,sizeof(std::size_t),1,pout);
      // chain list
      fwrite(list,sizeof(std::size_t),num,pout);
      // data
      const std::size_t dsize=6*num-3;
      double *dtemp=new double[dsize+1];
      dtemp[dsize]=dsize;
      backup(dtemp);
      fwrite(dtemp,sizeof(double),dsize+1,pout);
      // center-of-mass
      double cmass=cm.getMass();
      fwrite(&cmass,sizeof(double),1,pout);
      fwrite(cm.getPos(),sizeof(double),3,pout);
      fwrite(cm.getVel(),sizeof(double),3,pout);
      // mass of particles
      double *pmass=new double[num];
      p.getMassAll(pmass);
      fwrite(pmass,sizeof(double),num,pout);
      // interaction_parameters
      if (Int_pars) fwrite(Int_pars,sizeof(int_par),1,pout);

      fclose(pout);

      std::cerr<<"Chain dumping:\nNumber of stars ="<<num<<std::endl;
      std::cerr<<"Chain list index =";
      for (std::size_t i=0; i<num;i++) std::cerr<<list[i]<<" ";
      std::cerr<<"\nChain parameters t, B, w, X[][], V[][]= ";
      for (std::size_t i=0; i<dsize+1;i++) std::cerr<<dtemp[i]<<" ";
      std::cerr<<"\nMass of particles =";
      for (std::size_t i=0; i<num;i++) std::cerr<<pmass[i]<<" ";
      std::cerr<<std::endl;
      
      delete[] dtemp;
      delete[] pmass;
    }
  }

  //! Dump all data including chain, chainpars and necessary information (chain.smpars and #ds)
  /*!
   */
  

  //! Data reading from file function
  /*! Read chain data from binary file generated by dump(). 
      The first variable in the dumping data is particle number #num
      The second is the particle chain list index array #list
      The third is one dimensional array data storing #t, #B, #w, #X, #V generated by using backup().
      The forth is center of mass particle data
      The fifth is one dimensional array storing the masses of particles

      Notice the particle list will allocate new memory to store the particle data. The particle data are in the center-of-mass frame
      @param [in] filename: file for reading the data
   */
  void load(const char* filename) {
    std::FILE* pin = std::fopen(filename,"r");
    if (pin==NULL) std::cerr<<"Error: filename "<<filename<<" cannot be open!\n";
    else {
      // number of particles
      std::size_t n;
      std::size_t rn;
      rn = fread(&n,sizeof(std::size_t),1,pin);
      if(rn<1) {
        std::cerr<<"Error: cannot read number of particles!\n";
        abort();
      }
      if(n>nmax||p.getN()>0) {
        clear();
        allocate(n);
      }
      num = n;

      // chain list
      rn = fread(list,sizeof(std::size_t),n,pin);
      if(rn<n) {
        std::cerr<<"Error: reading chain list fails, required number of data is "<<n<<", only got "<<rn<<"!\n";
        abort();
      }
      
      // data
      const std::size_t dsize=6*n-3;
      double *dtemp=new double[dsize+1];
      rn = fread(dtemp,sizeof(double),dsize+1,pin);
      if(rn<=dsize) {
        std::cerr<<"Error: reading chain data fails, required data size is "<<dsize+1<<", only got "<<rn<<"!\n";
        abort();
      }
      
      restore(dtemp);

      // center-of-mass
      double cmass,cx[3],cv[3];
      rn = fread(&cmass,sizeof(double),1,pin);
      if(rn<1) {
        std::cerr<<"Error: cannot read center-of-mass!\n";
        abort();
      }
      rn = fread(cx,sizeof(double),3,pin);
      if(rn<3) {
        std::cerr<<"Error: reading position of center-of-mass particle fails, required reading 3 values, only got "<<rn<<"!\n";
        abort();
      }
      rn = fread(cv,sizeof(double),3,pin);
      if(rn<3) {
        std::cerr<<"Error: reading velocity of center-of-mass particle fails, required reading 3 values, only got "<<rn<<"!\n";
        abort();
      }

      cm.setMass(cmass);
      cm.setPos(cx[0],cx[1],cx[2]);
      cm.setVel(cv[0],cv[1],cv[2]);
      
      double *pmass=new double[n];
      rn = fread(pmass,sizeof(double),n,pin);
      if(rn<n) {
        std::cerr<<"Error: reading particle masses fails, required number of data is "<<n<<", only got "<<rn<<"!\n";
        abort();
      }

      Int_pars=new int_par;
      rn = fread(Int_pars,sizeof(int_par),1,pin);
      if (rn<1) {
        delete Int_pars;
        Int_pars=NULL;
      }
      
      F_load=true;
                 
      // mass of particles
      p.allocate(n);
      for (std::size_t i=0;i<n;i++) p[i].setMass(pmass[i]);

      // initialization
      resolve_XV();

      // set relative distance matrix, acceleration, potential and transformation parameter
      calc_rAPW();
      
      // kinetic energy
      calc_Ekin();

      F_Porigin = 0;
      F_Pmod = false;
      
      fclose(pin);

      std::cerr<<"Chain loading:\nNumber of stars ="<<n<<std::endl;
      std::cerr<<"Chain list index =";
      for (std::size_t i=0; i<n;i++) std::cerr<<list[i]<<" ";
      std::cerr<<"\nChain parameters t, B, w, X[][], V[][]= ";
      for (std::size_t i=0; i<dsize+1;i++) std::cerr<<dtemp[i]<<" ";
      std::cerr<<"\nMass of particles =";
      for (std::size_t i=0; i<n;i++) std::cerr<<pmass[i]<<" ";
      std::cerr<<std::endl;
      
      delete[] dtemp;
      delete[] pmass;

#ifdef TIME_PROFILE
      profile.initstep(pars->exp_itermax+1);
#endif
    }
  }

  //! Read particle data from file
  /*! Read particle data using fread (binary format). Number of particles should be first variable, then the data of particles (the data structure should be consistent with particle class).
      Notice if read is used, the particle data memory is allocated inside chainlist of chain.
    @param [in] filename: file to read the data
   */
  void readP(const char* filename) {
    p.read(filename);
  }
  
  //! Allocate memory for perturber list
  /*! Allocate memory for perturber particle list with maximum number of \a n
     @param [in] n: maximum number of perturbers
   */
  void initPext(const std::size_t n) {
    if (pext.getNmax()) {
      std::cerr<<"Error: Perturber list is already initialized!\n";
      abort();
    }
    pext.init(n);
  }

  
  //! Add one perturber particle
  /*! Add one perturber particle (address pointer) into #pext
     @param [in] a: perturber particle
   */
  void addPext(particle &a) { pext.add(a);}

  //! Add one chain as a perturber particle
  /*! Add one chain as a perturber particle (address pointer) into #pext
     @param [in] a: the chain perturber 
   */
  void addPext(chain<particle,int_par> &a) { pext.add(a);}

  //! Add a list of perturber particles
  /*! Add a list of perturber particles (address pointer) into #pext
     @param [in] n: number of particles need to be added
     @param [in] a: array of perturbers (address)
   */
  void addPext(const std::size_t n, particle a[]) { pext.add(n,a); }

  //! Remove one perturber 
  /*! Remove one perturber from #pext
     @param [in] i: perturber index in #pext needs to be removed
     @param [in] option: position update option: 
                 - true: shift last particle to current position (defaulted);
                 - false: shift all right particle to left by one
  */
  void removePext(const std::size_t i, bool option=true) { pext.remove(i,option); }
  
  //! Indicator of changes in particle list #p
  /*! Reture true if particle list is modifed, in this case, the chain may need to be initialized again
      \return true: Particle list is modified
  */
  bool isPmod() const { return F_Pmod; }

  //! Indicator of particle frame
  /*! Return true if particles positions and velocities in #p are in original frame
     \return: - 1: In the original frame
              - 2: Only positions are in the original frame
              - 0: In the center-of-mass frame
  */
  int isPorigin() const { return F_Porigin; }

  //! Get particle i (const reference)
  /*!
    \return: particle i (const reference)
   */
  const particle& getP(const std::size_t i) const {
    return p[i];
  }

  //! Get pertuber particle i (const reference)
  /*!
    \return: pertuber particle i (const reference)
   */
  const particle& getPext(const std::size_t i) const {
    return pext[i];
  }

  //! Get current number of chain members
  /*!
    Notice this is not always the number of particles.
    If AddP(), removeP() are used, init() should be called first to get consistent number.
    \return: number of chain members
  */
  const int getN() const{
    return num;
  }

  //! Initialization
  /*! Initialize chain based on particle list #p. After this function, positions and velocities of particles in #p will be shifted to their center-of-mass frame. \n
      Chain order list #list, relative position #X, velocity #V, initial system energy #Pt and initial time transformation parameter #w are calculated.
      The particle modification indicator (isPmod()) will be set to false.
    @param [in] time: current time of particle system
  */
  void init(const double time) {
    // update number indicator
    // update_num(p.getN());
    num = p.getN();
    
    // Generate chain link list
    generate_list();

    // set center-of-mass
    if (F_Porigin==1) {
      center_shift_init();
      F_Porigin=0;
    }
    else {
      std::cerr<<"Error: particles are not in original frame!\n";
      abort();
    }

    // set member relative position and velocity
    calc_XV();
    
    // set relative distance matrix, acceleration, potential and transformation parameter
    calc_rAPW();

    // Initial intgrt value t
    t = time;

    // kinetic energy
    calc_Ekin();

    // initial time step parameter
    Pt = -Pot - Ekin;
    w = W;

    // set F_Pmod to false
    F_Pmod = false;

#ifdef TIME_PROFILE
    profile.initstep(pars->exp_itermax+1);
#endif
    /*#ifdef DEBUG
//    for (int i=0;i<num;i++) {
//      std::cout<<i<<" m"<<std::setw(WIDTH)<<p[i].getMass()<<std::setw(WIDTH)<<"x";
//      for (int k=0;k<3;k++) std::cout<<std::setw(WIDTH)<<p[i].pos[k];
//      std::cout<<std::setw(WIDTH)<<"v ";
//      for (int k=0;k<3;k++) std::cout<<std::setw(WIDTH)<<p[i].vel[k];
//      std::cout<<std::endl;
//    }
//    std::cout<<std::setw(WIDTH)<<"kin:"<<std::setw(WIDTH)<<0.5*(p[0].getMass()*p[1].getMass()/(p[0].getMass()+p[1].getMass()))*(V[0][0]*V[0][0]+V[0][1]*V[0][1]+V[0][2]*V[0][2])<<std::endl;
    std::cerr<<std::setw(WIDTH)<<t;
    for (int i=0;i<num-1;i++) {
      for (int k=0;k<3;k++) std::cerr<<std::setw(WIDTH)<<X[i][k];
      for (int k=0;k<3;k++) std::cerr<<std::setw(WIDTH)<<V[i][k];
    }
    std::cerr<<std::setw(WIDTH)<<Ekin<<std::setw(WIDTH)<<Pot<<std::setw(WIDTH)<<Pt+Ekin-Pot<<std::setw(WIDTH)<<Pt<<std::setw(WIDTH)<<w<<std::endl;
    #endif*/
    
  }

  //! link interaction parameter class
  /*! link interaction parameter clsss to chain for the acceleration calculation functions
    @param[in] par: int_par type interaction parameter class. The address will be stored (not copy)
   */
  void link_int_par(int_par &par) {
    if (F_load) {
      std::cerr<<"Error: interaction parameter variable Int_pars is already set during load(), linking is forbidden!\n";
      abort();
    }
    else Int_pars = &par;
  }

  //! Inversed center-of-mass frame shift for particles
  /*! Shift the position and velocities of particles in #p from center-of-mass frame to original frame
      Notice the center-of-mass position and velocity use values from #cm
  */
  void center_shift_inverse() {
    if (F_Porigin==1) {
      std::cerr<<"Warning: particles are already in original frame!\n";
    }
    else {
      const double *rc = cm.getPos();
      const double *vc = cm.getVel();
      for (std::size_t i=0;i<num;i++) {
        if (F_Porigin==0) {
          const double *ri = p[i].getPos();
          p[i].setPos(ri[0] + rc[0],
                      ri[1] + rc[1],
                      ri[2] + rc[2]);
        }
        const double *vi = p[i].getVel();
        p[i].setVel(vi[0] + vc[0],
                    vi[1] + vc[1],
                    vi[2] + vc[2]);
      }
      F_Porigin = 1;
    }
  }

  //! Center-of-mass frame shift for particles
  /*! Shift positions and velocities of particles in #p from their original frame to their center-of-mass frame\n
      Notice the center-of-mass position and velocity use values from #cm
  */
  void center_shift() {
    if (F_Porigin>0) {
      const double *rc = cm.getPos();
      const double *vc = cm.getVel();
      for (std::size_t i=0;i<num;i++) {
        const double *ri = p[i].getPos();
        p[i].setPos(ri[0] - rc[0],
                    ri[1] - rc[1],
                    ri[2] - rc[2]);
        if (F_Porigin==1) {
          const double *vi = p[i].getVel();
          p[i].setVel(vi[0] - vc[0],
                      vi[1] - vc[1],
                      vi[2] - vc[2]);
        }
      }
      F_Porigin = 0;
    }
    else {
      std::cerr<<"Warning: particles are already in the center-of-mass frame!\n";
    }      
  }

  //! Backup chain data (#t, #Pt, #w, #X, #V)
  /*! Backup chain data to one dimensional array. 
      - #t : current physical time
      - #Pt : current time momentum (system binding energy)
      - #w : current time transformation parameter
      - #X : current relative position array
      - #V : current relative velocite array
     @param [out] db: backup array (size should be 6*#num-3) where #num is the total number of particles in #p
   */
  void backup(double* db) {
    const std::size_t dsize=6*num-3;
    if ((std::size_t)db[dsize]!=dsize) {
      std::cerr<<"Error: data array size ("<<(std::size_t)db[dsize]<<") for backup is not matched, should be ("<<dsize<<")!\n";
      abort();
    }
    const std::size_t ndata=3*(num-1);
    db[0] = t;
    db[1] = Pt;
    db[2] = w;
    std::memcpy(&db[3], X[0], ndata*sizeof(double));
    std::memcpy(&db[3+ndata], V[0], ndata*sizeof(double));
  }
     
  //! Restore chain data (#t, #Pt, #w, #X, #V)
  /*! Restore integration data from one dimensional array, the order of data should be #t, #Pt, #w, #X[#num][3], #V[#num][3]
      - #t : current physical time
      - #Pt : current time momentum (system binding energy)
      - #w : current time transformation parameter
      - #X : current relative position array
      - #V : current relative velocite array
     @param [in] db: one dimensional array that storing chain data (array size should be 6*#num-3) where #num is the total number of particles in #p
   */
  void restore(double* db) {
    const std::size_t dsize=6*num-3;
    if ((std::size_t)db[dsize]!=dsize) {
      std::cerr<<"Error: data array size ("<<(std::size_t)db[dsize]<<") for restore is not matched, should be ("<<dsize<<")!\n";
      abort();
    }
    const std::size_t ndata=3*(num-1);
    t = db[0];
    Pt = db[1];
    w = db[2];
    std::memcpy(X[0], &db[3], ndata*sizeof(double));
    std::memcpy(V[0], &db[3+ndata], ndata*sizeof(double));
  }

  //! Leapfrog integrator
  /*! Integration with Leapfrog method. \n
      The positions and velocities of particles in #p will be integrated in the center-of-mass frame
      @param [in] s: Integration step size
      @param [in] n: number of sub-steps needed to be divided. Integration step size is (\a s/\a n) and do \a n times as: X(s/2n)V(s/n)X(s/n)V(s/n)..X(s/2n)
      @param [in] force: external force (not perturber forces which are calculated in pert_force)
      @param [in] check_flag: 2: check link every step; 1: check link at then end; 0 :no check
      @param [in] dpoly: two dimensional array for storing 0 to \a (ndmax-*)'th order central difference of physical time #t at \a s/2 as a function of \a s, array size should be [ndmax][1]
      @param [in] ndmax: dpoly array size and the maximum difference is ndmax-*
  */             
//               recur_flag: flag to determine whether to resolve sub-chain particles for force calculations. notice this require the sub-chain to be integrated to current physical time. Thus is this a recursion call (tree-recusion-integration)
//               upforce: void (const particle * p, const particle *pext, double3* force). function to calculate force based on p and pext, return to force
  void Leapfrog_step_forward(const double s, const int n, const double3* force=NULL, int check_flag=1, double** dpoly=NULL, const int ndmax=0 ) {
#ifdef TIME_PROFILE
    profile.t_lf -= get_wtime();
#endif
    // Error check
    if (n<=0) {
      std::cerr<<"Error: step number shound be positive, current number is "<<n<<std::endl;
      abort();
    }
    /*    if (s<0) {
      std::cerr<<"Error: step size should be positive, current value is "<<s<<std::endl;
      abort();*/
    if (s==0) {
      std::cerr<<"Error: step size should not be zero!\n"<<std::endl;
      abort();
    }
    if (F_Porigin) {
      std::cerr<<"Error: particles are not in the center-of-mass frame, the integration can be dangerous!"<<std::endl;
      abort();
    }
    if (F_Pmod) {
      std::cerr<<"Error: particles are modified, initialization required!"<<std::endl;
      abort();
    }

    double ds = s/double(n);
    double3* ave_v=new double3[num];  // average velocity
    bool fpf = false;                 // perturber force indicator
    const int np = pext.getN();
    if (np>0) fpf = true;

    // for polynomial coefficient calculation first point
    // middle difference (first array is used to store the f_1/2)
    if(dpoly!=NULL) {
#ifdef TIME_PROFILE
      profile.t_dense -= get_wtime();
#endif
      if (pars->exp_sequence==3) {
        dpoly[1][0]=calc_dt_V(1.0);  //storage left edge derivate
        mid_diff_calc(&dpoly[4],ndmax-4,0,n);
      }
      else {
        dpoly[0][0]=calc_dt_V(1.0);  //storage left edge derivate
        edge_diff_calc(&dpoly[1],ndmax-1,0,n);
      }
#ifdef TIME_PROFILE
      profile.t_dense += get_wtime();
#endif
    }
                                               
    
    // integration loop
    for (std::size_t i=0;i<n;i++) {
      // half step forward for t (dependence: Ekin, Pt, w)
      double dt = calc_dt_X(ds*0.5);

      // step_forward time
      t += dt;
      
      // recursive integration-----------------------------------------------//
      /*
#ifdef TIME_PROFILE
      profile.t_lf += get_wtime();
#endif
      // check sub-chain 
      const std::size_t nct = p.getNchain();
      if (recur_flag&&nct>0) {
        // get chainlist table
        chain<particle,int_par> **clist = new chain<particle,int_par>*[nct];

        // looping to link sub-chains
        const std::size_t nt = p.getN();

        std::size_t k = 0;
        for (std::size_t j=0; j<nt; j++) {
          if (p.isChain(j)) {
            clist[k] = p.getSub(j);
            k++;
          }
        }

        // call tree recursion integration
        for (std::size_t j=0; j<k; j++) {
          // recursively call leapfrog until reaching the deepest branch.
          if (clist[j]->p.getNchain()>0)  clist[j]->Leapfrog_step_forward(ds, n, t, force, 1, true);
          else clist[j]->extrapolation_integration(ds, t, force);
        }
      }
#ifdef TIME_PROFILE
      profile.t_lf -= get_wtime();
#endif
      */
      //---------------------------------------------------------------------//

      // half step forward X (dependence: V, dt)
      step_forward_X(dt);

      // resolve X to p.v (dependence: X, cm.getMass())
      resolve_X();
      
      
      // perturber force
      if (fpf) {
        // get original position first, update p.x (dependence: X, cm.x)
        center_shift_inverse_X();
        // Update perturber force pf (dependence: pext, original-frame p.x, p.getMass())
        pert_force();
        // reset position to center-of-mass frame, update p.x 
        center_shift_X();
      }

      // Update rjk, A, Pot, dWdr, W for half X (dependence: pf, force, p.m, p.x, X)
      calc_rAPW(force); 

      // Update chain list order if necessary, update list, X, V (dependence: p.x, X, V)
      if (num>2&&check_flag==2) update_link();

      // Get time step dt(V) (dependence: Pot, W)
      double dvt = calc_dt_V(ds);

      // Step forward V (dependence: dt(V), V, A)
      step_forward_V(dvt);
      
      // Get averaged velocity, update p.x, p.v, ave_v (dependence: X, V)
      resolve_XV(ave_v);

      // forward Pt and w (dependence: dt(V), ave_v, p.m, p.v, force, dWdr, pf)
      step_forward_Ptw(dvt,ave_v,force,fpf);

      // Calcuale Kinetic energy (dependence: p.m, p.v)
      calc_Ekin();

      // step forward for t (dependence: Ekin, Pt, w)
      dt = calc_dt_X(ds*0.5);
      t += dt;

      // step forward for X (dependence: X, V)
      step_forward_X(dt);
      
      // for interpolation polynomial coefficient (difference)
      // middle difference (first array is used to store the f_1/2)
      if(dpoly!=NULL) {
#ifdef TIME_PROFILE
        profile.t_dense -= get_wtime();
#endif
        if(pars->exp_sequence==3) {
          if (i==n/2-1) {
            dpoly[0][0]=t; // y
            dpoly[3][0]=2.0*dt/ds; // f(x)
//#ifdef DEBUG
//            std::cerr<<"Mid time = "<<t<<", n="<<n<<"; i="<<i+1<<std::endl;
//#endif
          }
          mid_diff_calc(&dpoly[4],ndmax-4,i+1,n);
        }
        // edge difference
        else edge_diff_calc(&dpoly[1],ndmax-1,i+1,n);
#ifdef TIME_PROFILE
        profile.t_dense += get_wtime();
#endif
      }
    }

    // resolve X at last, update p.x (dependence: X)
    resolve_X();

    // Update rjk, A, Pot, dWdr, W (notice A will be incorrect since pf is not updated)
    calc_rAPW(force);

#ifdef TIME_PROFILE
    profile.t_dense -= get_wtime();
#endif
    if(dpoly!=NULL) {
      if(pars->exp_sequence==3) dpoly[2][0]=calc_dt_V(1.0);
      else dpoly[0][1]=calc_dt_V(1.0);
    }
#ifdef TIME_PROFILE
    profile.t_dense += get_wtime();
#endif

//#ifdef DEBUG
//    std::cerr<<"Ending time = "<<t<<", n="<<n<<std::endl;
//#endif
//#ifdef DEBUG
//    std::cerr<<std::setw(WIDTH)<<t;
//    //      for (int i=0;i<num;i++) 
//    //        for (int k=0;k<3;k++) std::cerr<<p[i].pos[k]<<" ";
//    for (int i=0;i<num-1;i++) {
//      for (int k=0;k<3;k++) std::cerr<<std::setw(WIDTH)<<X[i][k];
//      for (int k=0;k<3;k++) std::cerr<<std::setw(WIDTH)<<V[i][k];
//    }
//    std::cerr<<std::setw(WIDTH)<<Ekin<<std::setw(WIDTH)<<Pot<<std::setw(WIDTH)<<Pt+Ekin-Pot<<std::setw(WIDTH)<<Pt<<std::setw(WIDTH)<<w<<std::endl;
//#endif
      
    // update chain list order and calculate potential
    if(num>2&&check_flag==1)  update_link();

    // clear memory
    delete[] ave_v;
#ifdef TIME_PROFILE
    profile.t_lf += get_wtime();
#endif
  }

  //! Extrapolation integration
  /*! Use extrapolation method to get highly accurate integration based on Leapfrog_step_forward().
    The auto-determination of extrapolation orders based on the accuracy requirement is used. 

     @param [in] ds: integration step size
     @param [in] toff: ending physical time
                      - if value is negative, it means integration will be done with fixed step size \a ds
                      - if value is positive and after step \a ds, the ending physical time is larger than \a toff, the interpolation of physical time #t (dense output) will be done instead of integration. In this case, the data are kept as initial values. Instead, the returning value is the ds modification factor (negative value), which can be used to modified current \a ds and redo the integration by calling this function again with new ds to approach ending physical time of \a toff. Notice if the required time sychronization criterion (set in chainpars.setEXP()) is small (<phase and energy error criterion), several iteration may be needed to get the physical time below this criterion.
     @param [in] force: external force (not perturber forces which are calculated in pert_force)
     @param [in] err_ignore: if true, force integration and ignore error criterion (default false)
     \return factor
            - if factor is positive, it is optimized step size modification factor for next step (\a ds *= factor)
            - if factor is negative and \a toff>0; the -factor is used for calculate new ds' = -factor * \a ds. Thus this function should be called again with new step size ds' and the new result should have ending physical time close to \a toff.
            - if factor is zero, maximum extrapolation sequence index (accuracy order/iteration times) is fixed (\ref ARC::chainpars) and err_ingore is false, it means the error criterion cannot be satisfied with current maximum sequence index. In this case no integration is done and the data are kept as initial values. User should reduce the integration step and re-call this function.
   */
  double extrapolation_integration(const double ds, const double toff=-1.0, const double3* force=NULL, const bool err_ignore=false) {
#ifdef TIME_PROFILE
    profile.t_ep -= get_wtime();
#endif
    // get parameters
    const double error = pars->exp_error;
    const std::size_t itermax = pars->exp_itermax;
    const int method = pars->exp_method;
    const int sq = pars->exp_sequence;
    const int *step = pars->step;
    
    // array size indicator for relative position and velocity
    const std::size_t nrel = num-1;
    // data storage size (extra one is used to show the size of the array for safety check)
    const std::size_t dsize = 6*nrel+3;
    // edge difference array size;
    // const std::size_t psize=dsize*2;
    
    // for storage
    // for convenient, the data are storaged in one array with size (2*nrel+1)*3, which represents t, Pt, w, X[3][nrel], V[3][nrel]
    double d0[dsize+1],dtemp[dsize+1];
    double* dn[itermax];
    d0[dsize] = (double)dsize;    // label for safety check
    dtemp[dsize] = (double)dsize; // label for safety check
    for (std::size_t i=0; i<itermax; i++) {
      dn[i] = new double[dsize+1];
      dn[i][dsize] = (double)dsize; // label for safety check
    }
    double Ekin0,Pot0;

    // for dense output polynomial
#ifdef TIME_PROFILE
    profile.t_dense -= get_wtime();
#endif
    bool ip_flag = true;  // interpolation coefficient calculation flag
    double** pd[itermax]; // central difference, [*] indicate different accuracy level 
    int ndmax[itermax];   // maximum difference order
    std::size_t pnn;      // data size
    if(sq==3) pnn = 1;     // middle difference case
    else pnn = 2;         // edge two points case
#ifdef TIME_PROFILE
    profile.t_dense += get_wtime();
#endif

    // for error check
    double3 CX={};
    double3 CXN;
    double cxerr=error+1.0;
    double eerr=error+1.0;
    double cxerr0=cxerr+1.0;
    double eerr0=eerr+1.0;
    // error for step estimation
    double werrmax=std::numeric_limits<double>::max();

    // backup initial values (t, Pt, w, X, V, Ekin)
    backup(d0);
    Ekin0 = Ekin;
    Pot0 = Pot;

    // new step
    double dsn = 1.0;
    
    std::size_t intcount = 0;  // iteration counter
    std::size_t itercount = 0; // iteration efforts count
    
    
    while (true) {
      // convergency check
      if (!pars->exp_fix_iter) {
        if (cxerr <= 0.5*error && std::abs(eerr) <= 0.1*error) break;
        
        if (cxerr>=cxerr0 || std::abs(eerr) >= std::abs(eerr0)) {
          if (cxerr < error && std::abs(eerr) < error) break;
          else if (intcount > std::min((size_t)10,itermax)){
#ifdef ARC_WARN            
            std::cerr<<"Warning: extrapolation cannot converge anymore, energy error - current: "<<eerr<<"  previous: "<<eerr0<<"   , phase error - current: "<<cxerr<<"  previous: "<<cxerr0<<", try to change the error criterion (notice energy error is cumulative value)\n";
#endif
            // in the case of serious energy error, quit the simulation and dump the data
            if (std::abs(eerr)*std::min(1.0,std::abs(Pt))>100.0*error) {
              std::cerr<<"Error!: extrapolation cannot converge anymore, but energy error is too large. energy error - current: "<<eerr<<"  previous: "<<eerr0<<"   , phase error - current: "<<cxerr<<"  previous: "<<cxerr0<<", ds = "<<ds<<", N="<<num<<", the particle data before integration is dumped to file \"ARC_dump.dat\"\n";
              restore(d0);
              dump("ARC_dump.dat");
              pars->dump("ARC_dump.par");
              
              abort();
            }
            break;
          }
        }
        // if completely converged, check energy error
        if (cxerr==0) {
          if (std::abs(eerr) > error)
#ifdef ARC_WARN            
            std::cerr<<"Warning: phase error reach zero but energy error "<<eerr<<" cannot reach criterion "<<error<<"!\n";
#endif
          // in the case of serious energy error, quit the simulation and dump the data
          if (std::abs(eerr)>1000.0*error) {
            std::cerr<<"Error!: extrapolation cannot converge anymore, but energy error is too large. energy error - current: "<<eerr<<"  previous: "<<eerr0<<"   , phase error - current: "<<cxerr<<"  previous: "<<cxerr0<<", the particle data before integration is dumped to file \"ARC_dump.dat\"\n";
            restore(d0);
            dump("ARC_dump.dat");
            pars->dump("ARC_dump.par");
            abort();
          }
          break;
        }
      }

      // iteration limit check
      if (intcount == itermax) {
        if(cxerr < error && std::abs(eerr) < error) break;
        else if(err_ignore) break;
        else if(pars->exp_fix_iter) {
          // reset the initial data (t, Pt, w, X, V, Ekin)
          restore(d0);
          Ekin = Ekin0;
          Pot  = Pot0;
          // reset velocity to get correct w
          if (pars->beta>0) resolve_V();

          // clear memory
          for (std::size_t i=0; i<itermax; i++) {
            delete[] dn[i];
          }
          for (std::size_t i=0; i<intcount-1; i++) {
            if (pd[i]!=NULL) {
              for (std::size_t j=0; j<ndmax[i]; j++) delete[] pd[i][j];
              delete[] pd[i];
            }
          }

#ifdef TIME_PROFILE
          profile.stepcount[intcount]++;
          profile.t_ep += get_wtime();
          profile.itercount +=itercount;
#endif
          // if error is big, reduce the integration step and recursive call extrapolation
          //double dsfactor = std::min(EP::H_opt_factor(std::max(cxerr,std::abs(eerr)),error,intcount),0.01);
          //std::cerr<<"Accuracy not reached, time="<<t<<"; phase error="<<cxerr<<"; energy error="<<eerr<<"; error criterion="<<error<<"; try reducing step size to "<<ds*dsfactor<<std::endl;
          //extrapolation_integration(ds*dsfactor,toff,force);

          // indicate of error too large
          return 0.0;
        }
        if (std::abs(eerr) > error) {
          std::cerr<<"Error: maximum iteration step number "<<itermax<<" reached, but energy error "<<eerr<<" is larger than criterion "<<error<<std::endl;
        } else {
          std::cerr<<"Error: maximum iteration step number "<<itermax<<" reached, but phase error "<<cxerr<<" is larger than criterion "<<error<<std::endl;
        }          
        abort();
      }
      
      if (intcount>0) {
        // reset the initial data (t, Pt, w, X, V, Ekin)
        restore(d0);
        Ekin = Ekin0;
        Pot = Pot0;
        // reset velocity to get correct w
        if (pars->beta>0) resolve_V();
      }

      // Dense output
#ifdef TIME_PROFILE
      profile.t_dense -= get_wtime();
#endif
      if (ip_flag) {
        // middle difference case: difference order from 1 to 2*intcount+2 (2*kappa-2; kappa=intcount+1), first one is used to storage f(x), 2,3 are used for edge f'(x) at (0, ds)
        if(sq==3) ndmax[intcount] = 2*intcount+5;
        // edge difference case: difference order from 1 to intcount+1, first store f(x)
        else ndmax[intcount] = intcount+2;
        
        // pd[][*]: * storage f(x) and difference
        pd[intcount] = new double*[ndmax[intcount]];
        // pd[][][*]: * storage data (t, Pt, w, X, V)
        for (std::size_t j=0;j<ndmax[intcount];j++) pd[intcount][j] = new double[pnn];
      }
      else pd[intcount] = NULL;
#ifdef TIME_PROFILE
      profile.t_dense += get_wtime();
#endif

      // intergration
      Leapfrog_step_forward(ds,step[intcount],force,0,pd[intcount],ndmax[intcount]);

      // increase iteration counter
      itercount += step[intcount];
      
      if (intcount == 0) {
        // storage the results to [n]
        backup(dn[intcount]);

        // relative position vector between first and last particle for phase error check
        for (std::size_t i=0;i<3;i++) CX[i] = 0.0;
        for (std::size_t i=0;i<num-1;i++) {
          CX[0] += X[i][0];
          CX[1] += X[i][1];
          CX[2] += X[i][2];
        }
      }
      else {
        // storage the results to [temp]
        backup(dtemp);
        
        // iteration
        // Using Polynomial method
        if (method==1) EP::polynomial_extrapolation(dn,dtemp,step,dsize,intcount);
        // Using Rational interpolation method
        else EP::rational_extrapolation(dn,dtemp,step,dsize,intcount);

        // set final results back to chain array
        restore(dn[intcount]);
 
        // resolve particle
        resolve_XV();
        // recalculate the energy
        calc_Ekin();
        // force, potential and W
        calc_rAPW(force);

        // phase error calculation
        for (std::size_t i=0;i<3;i++) CXN[i] = 0.0;
        for (std::size_t i=0;i<nrel;i++) {
          CXN[0] += X[i][0];
          CXN[1] += X[i][1];
          CXN[2] += X[i][2];
        }
        double RCXN2 = CXN[0]*CXN[0] + CXN[1]*CXN[1] + CXN[2]*CXN[2];
      
        double dcx1 = CXN[0] - CX[0];
        double dcx2 = CXN[1] - CX[1];
        double dcx3 = CXN[2] - CX[2];
        cxerr0 = cxerr;
        cxerr = std::sqrt((dcx1*dcx1 + dcx2*dcx2 + dcx3*dcx3)/RCXN2);
        eerr0 = eerr;
        eerr = (Ekin+Pot+Pt-Ekin0-Pot0-d0[1])/Pt;
        //        std::cerr<<"Ekin="<<Ekin<<" Pot="<<Pot<<" Pt="<<Pt<<" Ekin0="<<Ekin0<<" Pot0="<<Pot0<<" Pt0="<<d0[1]<<" eerr="<<eerr<<std::endl;
        std::memcpy(CX,CXN,3*sizeof(double));

        if (pars->auto_step==1) {
          // get error estimation
          double ermax=std::min(EP::extrapolation_error(dn,dsize,intcount),std::min(eerr,cxerr));
          double dsfactor = EP::H_opt_factor(ermax,error,intcount+1);
          double werrn = ((double)itercount+num)/ dsfactor;
          if (ermax>0&&werrn<werrmax) {
            werrmax = werrn;
            dsn = dsfactor;
#ifdef DEBUG
            std::cerr<<"ERR factor update: sequence="<<step[intcount]<<"; modify factor="<<dsfactor<<"; ermax="<<ermax<<"; eerr="<<eerr<<"; cxerr="<<cxerr<<"; ds="<<ds<<std::endl;
#endif
          }
        }
      }
      
      intcount++;
#ifdef DEBUG
      std::cerr<<std::setprecision(6)<<"Iter.= "<<intcount<<" Dep.= "<<step[intcount]<<" P-err.= "<<cxerr;
      std::cerr<<" E-err.="<<Ekin+Pot+Pt-Ekin0-Pot0-d0[1]<<" Pt ="<<std::setprecision(12)<<Pt<<std::endl;
#endif 
    }

    // for dense output
    if (ip_flag&&toff>0&&toff<t&&std::abs(toff-t)>pars->dterr) {
#ifdef TIME_PROFILE
      profile.t_dense -= get_wtime();
#endif

#ifdef DEBUG
      std::cerr<<"ds="<<ds<<" step[0]="<<step[0]<<" terr="<<toff-t<<" t="<<t<<" toff="<<toff<<std::endl;
#endif
      // calculate derivates from differences
      for (std::size_t i=0; i<intcount; i++) {
        // first difference pointer
        double **pdptr;
        int dpsize;
        double h;
        if(sq==3) {
          // middle difference case first element is f(x)
          pdptr=&pd[i][4];
          // differece order number should be reduced by one
          dpsize = ndmax[i]-4;
          // step size
          h = 2*ds/(double)step[i];
        }
        else {
          pdptr=&pd[i][1];
          dpsize = ndmax[i]-1;
          h = ds/(double)step[i];
        }

        diff_dev_calc(pdptr,h,dpsize,pnn); 
      }        
      
      // extrapolation table
      double* pn[intcount];
      for (std::size_t i=0; i<intcount; i++) pn[i] = new double[pnn];

      // starting accuracy index
      std::size_t istart=0;
      
      // from low difference order (1) to high difference order (ndmax) 
      for (std::size_t i=0; i<(std::size_t)ndmax[intcount-1]; i++) {

        // find correct istart;
        if (i>=ndmax[istart]) istart++;

        if (istart>=intcount) break;
        
        // storage result of first order accuracy
        std::memcpy(pn[0],pd[istart][i],pnn*sizeof(double));
#ifdef DEBUG
        std::cerr<<"Poly calc order="<<istart<<" step("<<istart<<")="<<step[istart]<<" t X11^("<<i+1<<")_"<<istart<<"="<<pd[istart][i][0]<<"\t"<<pd[istart][i][3]<<std::endl;
#endif
        // extrapolation to higher order accuracy
        for (std::size_t j=istart+1; j<intcount; j++) {
#ifdef DEBUG
          std::cerr<<"Poly calc order="<<j<<" step("<<istart<<")="<<step[istart]<<" t X11^("<<i+1<<")_"<<j<<"="<<pd[j][i][0]<<"\t"<<pd[j][i][3]<<std::endl;
#endif
          if (method==1) EP::polynomial_extrapolation(pn,pd[j][i],&step[istart],pnn,j-istart);
          if (method==2) EP::rational_extrapolation(pn,pd[j][i],&step[istart],pnn,j-istart);
#ifdef DEBUG
          std::cerr<<"Poly extra order="<<j-istart<<" step("<<istart<<")="<<step[istart]<<" t X11^("<<i+1<<")_"<<j<<"="<<pd[j][i][0]<<"\t"<<pd[j][i][3]<<std::endl;
#endif
        }
#ifdef DEBUG
          std::cerr<<"Final result t X11^("<<i+1<<")="<<pd[intcount-1][i][0]<<"\t"<<pd[intcount-1][i][3]<<std::endl;
#endif
      }

      // number of points
      int npoints;
      // store position s
      double* xpoint;
      // maximum difference level
      int* nlev;
      // polynomial
      double* pcoff;
      // data point
      double** fpoint;
      // final derivate starting pointer
      double*** dfptr;
      

      // for middle difference case (1 point)
      if(sq==3) {

        npoints=3;
        
        xpoint=new double[npoints];
        xpoint[0]=0.0;
        xpoint[1]=0.5*ds;
        xpoint[2]=ds;

        nlev=new int[npoints];
        nlev[0] = 2;
        nlev[1] = ndmax[intcount-1]-2;
        nlev[2] = 2;

        fpoint=new double*[npoints];
        fpoint[0] = d0;
        fpoint[1] = pd[intcount-1][0];
        fpoint[2] = dn[intcount-1];

        // \sum nlev = 2*intcount+2;
        pcoff = new double[ndmax[intcount-1]+2];

        dfptr=new double**[nlev[1]-1];
        for (int i=0;i<nlev[1]-1;i++) {
          dfptr[i]=new double*[3];
          dfptr[i][0]=NULL;
          dfptr[i][1]=pd[intcount-1][i+3];
          dfptr[i][2]=NULL;
        }
        dfptr[0][0]=pd[intcount-1][1]; //left edge 
        dfptr[0][2]=pd[intcount-1][2]; //right edge
      }
      // for edge difference case (2 points)
      else {
        npoints=2;

        xpoint=new double[npoints];
        xpoint[0]=0.0;
        xpoint[1]=ds;

        nlev=new int[npoints];
        nlev[0] = ndmax[intcount-1]+1;
        nlev[1] = nlev[0];
        
        // store f(x)
        fpoint=new double*[npoints];
        fpoint[0] = d0;
        fpoint[1] = dn[intcount-1];

        // \sum nlev = 2*intcount+2;
        pcoff= new double[2*ndmax[intcount-1]+2];

        dfptr=new double**[nlev[0]-1];
        for (int i=0;i<nlev[0]-1;i++) {
          dfptr[i] = new double*[2];
          dfptr[i][0] = pd[intcount-1][i];
          dfptr[i][1] = &pd[intcount-1][i][1];
        }
      }

      
      // Hermite interpolation
      EP::Hermite_interpolation_coefficients(&pcoff,xpoint,fpoint,dfptr,1,npoints,nlev);

#ifdef DEBUG
      std::cerr<<"PCOFF: ";
      for (std::size_t i=0;i<ndmax[intcount-1]+2;i++) std::cerr<<" "<<pcoff[i];
      std::cerr<<std::endl;
#endif

      // Iteration to get correct physical time position
      double dsi[2]   = {0,ds};    // edges for iteration
      double tsi[2]   = {d0[0],t}; // edges values
      double dsm,tpre;             // expected ds and t;
      const double dterr = pars->dterr;
      const double dterr3 = 1000*dterr;  // 1000 * time error criterion

      bool rf_method=false;
      int find_root_count=0;
      do {
        if (rf_method) {
          dsm = (dsi[0]*(tsi[1]-toff)-dsi[1]*(tsi[0]-toff))/(tsi[1]-tsi[0]);  // Use regula falsi method to find accurate ds
        }
        else {
          dsm = (dsi[0]+dsi[1])*0.5;      // Use bisection method to get approximate region
        }
        
        EP::Hermite_interpolation_polynomial(dsm,&tpre,&pcoff,xpoint,1,npoints,nlev);
        if (tpre > toff) {
          if (dsi[1]==dsm) break;
          dsi[1] = dsm;
          tsi[1] = tpre;
        }
        else {
          if (dsi[0]==dsm) break;
          dsi[0] = dsm;
          tsi[0] = tpre;
        }
#ifdef DEBUG
        std::cerr<<std::setprecision(15)<<"Find root: dsm="<<dsm<<"; t="<<tpre<<"; error="<<tpre-toff<<"; ds="<<ds<<std::endl;
#endif
        if(std::abs(tpre-toff)<dterr3) rf_method=true;
        find_root_count++;
        if (find_root_count>100) {
          //          std::cerr<<"Error! can not find ds to reach physical time toff. current searching ds="<<dsm<<"; current interpolated time="<<tpre<<"; toff="<<toff<<"; ds="<<ds<<"; err="<<tpre-toff<<std::endl;
          //          abort();
          break;
        }
      } while (std::abs(tpre-toff)>0.1*dterr);

      // Get the interpolation result
      //EP::Hermite_interpolation_polynomial(dsm,dtemp,&pcoff,xpoint,1,npoints,nlev);

      // update time factor
      dsn = -(dsm/ds);

      // avoid energy issue
      //if(cxerr < error && std::abs(eerr) < error) dsn = 0.5*ds;
      
      // update the results
      //restore(dtemp);

      // reset the data to original
      restore(d0);
      Ekin = Ekin0;
      Pot  = Pot0;
      // reset velocity to get correct w
      if (pars->beta>0) resolve_V();

#ifdef DEBUG
      /*std::cerr<<"Getting: ds= "<<dsm;
      for (std::size_t i=0;i<dsize;i++) std::cerr<<" "<<dtemp[i];
      std::cerr<<std::endl;
      
      EP::Hermite_interpolation_polynomial(0,dtemp,&pcoff,xpoint,dsize,2,nlev);
      std::cerr<<"Starting:";
      for (std::size_t i=0;i<dsize;i++) std::cerr<<" "<<dtemp[i];
      std::cerr<<std::endl;

      std::cerr<<"Ending:";
      EP::Hermite_interpolation_polynomial(ds,dtemp,&pcoff,xpoint,dsize,2,nlev);
      for (std::size_t i=0;i<dsize;i++) std::cerr<<" "<<dtemp[i];
      std::cerr<<std::endl;
      */
      for (std::size_t i=0; i<=1000; i++) {
        dsm = ds/1000*i;
        EP::Hermite_interpolation_polynomial(dsm,&tpre,&pcoff,xpoint,1,npoints,nlev);
//        EP::Hermite_interpolation_polynomial(dsm,dtemp,pcoff,xpoint,dsize,npoints,nlev);
//        // update the results
//        restore(dtemp);
//        // resolve particle
//        resolve_XV();
//        // recalculate the energy
//        calc_Ekin();
//        // force, potential and W
//        calc_rAPW(force);
//      
        std::cerr<<"Loop: "<<dsm;
        std::cerr<<" "<<std::setprecision(15)<<tpre;
//        std::cerr<<" "<<(Ekin-Pot+Pt)/Pt;
//        for (std::size_t i=0;i<dsize;i++) std::cerr<<std::setprecision(10)<<" "<<dtemp[i];
        std::cerr<<std::endl;
      }
      
#endif

      // clear memory
      for (std::size_t i=0; i<intcount; i++) delete[] pn[i];
      for (std::size_t i=0; i<(std::size_t) nlev[1]-1; i++) delete dfptr[i];
      delete[] dfptr;
      delete[] pcoff;
      delete[] xpoint;
      delete[] nlev;
      delete[] fpoint;

//      // clear memory
//      for (std::size_t i=0; i<itermax; i++) {
//        delete[] dn[i];
//      }
//      for (std::size_t i=0; i<intcount; i++) {
//        if (pd[i]!=NULL) {
//          for (std::size_t j=0; j<ndmax[i]; j++) delete[] pd[i][j];
//          delete[] pd[i];
//        }
//      }
//#ifdef TIME_PROFILE
//      profile.stepcount[intcount+1]++;
//      profile.t_ep += get_wtime();
//      profile.itercount +=itercount;
//#endif
//      // recursive calling extrapolation_integrator to get correct time
//      extrapolation_integration(dsm,toff,force);
//      
//      return dsn;
#ifdef TIME_PROFILE
      profile.t_dense += get_wtime();
#endif
    }
    else {
      // auto-step
      if      (pars->auto_step==1)
        dsn = std::min(std::max(dsn,pars->auto_step_fac_min),pars->auto_step_fac_max);
      else if (pars->auto_step==2) {
        dsn = calc_next_step_XVA()/ds;
        dsn = std::min(std::max(dsn,pars->auto_step_fac_min),pars->auto_step_fac_max);
      }
      else if (pars->auto_step==3) {
        if      (intcount>pars->auto_step_iter_max) dsn = pars->auto_step_fac_min;
        else if (intcount<pars->auto_step_iter_min) dsn = pars->auto_step_fac_max;
        else    dsn = 1.0;
      }
      else if (pars->auto_step==4)
        dsn = calc_next_step_custom()/ds;

      // update chain link order
      if(num>2) update_link();
    }
    
    // clear memory
    for (std::size_t i=0; i<itermax; i++) {
      delete[] dn[i];
    }
#ifdef TIME_PROFILE
    profile.t_dense -= get_wtime();
#endif
    if (ip_flag) {
      for (std::size_t i=0; i<intcount; i++) {
        if (pd[i]!=NULL) {
          for (std::size_t j=0; j<ndmax[i]; j++) delete[] pd[i][j];
          delete[] pd[i];
        }
      }
    }
#ifdef TIME_PROFILE
    profile.t_dense += get_wtime();
#endif

#ifdef TIME_PROFILE
    profile.stepcount[intcount+1]++;
    profile.t_ep += get_wtime();
    profile.itercount +=itercount;
#endif

    return dsn;
  }

  //! Get current physical time
  /*! \return current physical time
   */
  double getTime() const {
    return t;
  }
  //! Get current kinetic energy
  /*! \return current kinetic energy
  */
  double getEkin() const {
    return Ekin;
  }

  //! Get current potential energy
  /*! \return current potetnial energy (negative value for bounded systems)
  */
  double getPot() const {
    return Pot;
  }

  //! Get current time momemtum \f$Pt\f$ (current system binding energy)
  /*! \return time momemtum \f$Pt\f$ (current system binding energy \f$(-H(t))\f$)
  */
  double getPt() const {
    return Pt;
  }
  
  //! Get current integrated time transformation function value \f$w\f$
  /*! \f$w\f$: \f$ \frac{dw}{dt} = \sum_k \frac{\partial W}{\partial \vec{r_k}} \bullet \vec{v_k} \f$
      (see W in getW();)
      \return \f$w\f$
  */
  double getw() const {
    return w;
  }

  //! Get current time transformation function value \f$W\f$
  /*! The \f$W\f$ is defined in ::ARC::pair_AW().
      \return \f$W\f$
  */
  double getW() const {
    return W;
  }
  
#ifdef DEBUG
  /* set_X
     function: modify X
     argument: i: the index of X
               k: axis-0:x,1:y,2:z
               value: new value
  void set_X(const std::size_t i, const std::size_t k, const double value) {
    X[i][k] = value;
  }
   */


  //! [DEBUG] print chain data
  /*! Print chain data to std::cerr for debuging purpose 
      Only avaiable when macro name DEBUG is defined
     @param [in] width: digital width of output values
  */
  void print(const int width=15) {
    if (width<=0) {
      std::cerr<<"Error: width should be larger than zero!\n";
      abort();
    }
    char xyz[4]={'x','y','z','r'};
    std::cerr<<"---- particle list------\n ";
    for (std::size_t i=0;i<num;i++) std::cerr<<std::setw(width)<<list[i];
    std::cerr<<"\n----- relative position X ------\n";
    for (std::size_t k=0;k<3;k++) {
      std::cerr<<xyz[k];
      for (std::size_t i=0;i<num-1;i++) std::cerr<<std::setw(width)<<X[i][k];
      std::cerr<<std::endl;
    }
    std::cerr<<"\n----- relative velocity V ------\n";
    for (std::size_t k=0;k<3;k++) {
      std::cerr<<xyz[k];
      for (std::size_t i=0;i<num-1;i++) std::cerr<<std::setw(width)<<V[i][k];
      std::cerr<<std::endl;
    }
    std::cerr<<"\n----- Acceleration A ------\n";
    for (std::size_t k=0;k<3;k++) {
      std::cerr<<xyz[k];
      for (std::size_t i=0;i<num;i++) std::cerr<<std::setw(width)<<acc[i][k];
      std::cerr<<std::endl;
    }
    std::cerr<<"\n----- part omega / part rk ------\n";
    for (std::size_t k=0;k<3;k++) {
      std::cerr<<xyz[k];
      for (std::size_t i=0;i<num;i++) std::cerr<<std::setw(width)<<dWdr[i][k];
      std::cerr<<std::endl;
    }
    std::cerr<<"\n----- system parameters ------\n"
             <<"---Center-of-mass data: \n"
             <<"mass:"<<std::setw(width)<<cm.getMass()<<std::endl
             <<"(particle*)cm.pos:"
             <<std::setw(width)<<cm.getPos()[0]
             <<std::setw(width)<<cm.getPos()[1]
             <<std::setw(width)<<cm.getPos()[2]<<std::endl
             <<"(particle*)cm.vel:"
             <<std::setw(width)<<cm.getVel()[0]
             <<std::setw(width)<<cm.getVel()[1]
             <<std::setw(width)<<cm.getVel()[2]<<std::endl
             <<"---Energy: "<<std::endl
             <<"Kinetic energy Ekin:"<<std::setw(width)<<Ekin<<std::endl
             <<"Potential energy Pot:"<<std::setw(width)<<Pot<<std::endl
             <<"Transformation factor Omega:"<<std::setw(width)<<W<<std::endl
             <<"---Integration parameters: "<<std::endl
             <<"Time momentum Pt:"<<std::setw(width)<<Pt<<std::endl
             <<"Transformation coefficient omega:"<<std::setw(width)<<w<<std::endl
             <<"---Time step coefficients:"<<std::endl
             <<"alpha: "<<std::setw(width)<<pars->alpha<<std::endl
             <<"beta: "<<std::setw(width)<<pars->beta<<std::endl
             <<"gamma: "<<std::setw(width)<<pars->gamma<<std::endl;
  }
#endif  
};

//! Generalized list to store chain and particle members
/*! A list that storing chain and particle memory addresses (based on template class particle)
 */
template <class particle, class int_par>
class chainlist{
  std::size_t num; //!< number of current particles in the list #p
  std::size_t nmax; //!< maximum number of particles allown to store
  std::size_t nchain;  //!< number of chain type members
  bool* cflag; //!< flag array to indicate whether the corresponding particle in #p is chain (true: chain; false: Particle)
  void** p; //!< particle list array (void pointer array)
  bool alloc_flag ;//!< indicate whether the p stores the particle memory address (false) or allocates new memory (true)

public:
  //! Constructor 
  /*! Set current particle number to zero, need to use init() to allocate memory for storing particle addresses
   */
  chainlist(): num(0), nmax(0), nchain(0), alloc_flag(false) {};
  
  //! Constructor with maximum number of particle \a n
  /*! Set maximum particle number to \a n and allocate memory for #p to storing the particle addresses (maximum \a n)
   */
  chainlist(const std::size_t n): alloc_flag(false) { init(n); }

  //! Initialization with memory allocation
  /*! Set maximum particle number to \a n and allocate memory for #p to storing the particle addresses (maximum \a n)
   */
  void init(const std::size_t n) {
    num = 0;
    nmax = n;
    nchain = 0;
    cflag=new bool[n];
    p=new void*[n];
  }

  //! allocate memory of particles
  /*! allocate memory of n particle data with type of particle class
    @param [in] n: number of particles
   */
  void allocate(const std::size_t n) {
    if (num>0&&!alloc_flag) {
      std::cerr<<"Error: the particle list already stores particle addresses, the allocation of new particle memory is forbidden, please clear first!\n";
      abort();
    }
    if(n+num>nmax) {
      std::cerr<<"Error: particle number allocated ("<<n<<") + current particle number ("<<num<<") > maximum number allown ("<<nmax<<")!\n";
      abort();
    }
    for (std::size_t i=0; i<n; i++) {
      p[i+num] = new particle;
      cflag[i+num] = false;
    }
    num += n;
    alloc_flag = true;
  }
  
  //! Clear function
  /*! Free dynamical memory space used in particle address list #p
   */
  void clear() {
    if (nmax>0) {
      if (alloc_flag) 
        for (std::size_t i=0;i<num;i++)
          if (p[i]!=NULL) delete (particle*)p[i];
      nmax = 0;
      num = 0;
      delete[] cflag;
      delete[] p;
      alloc_flag = false;
    }
  }

  //! Destructor
  ~chainlist() {
    if (nmax>0) {
      if (alloc_flag) 
        for (std::size_t i=0;i<num;i++)
          if (p[i]!=NULL) delete (particle*)p[i];
      delete[] cflag;
      delete[] p;
    }
  }

  //! Get current particle number
  /*! \return Current particle number in particle address list #p
   */
  std::size_t getN() const {
    return num;
  }

  //! Get maximum particle number
  /*! \return Current maximum particle number that can be stored in particle address list #p
   */
  std::size_t getNmax() const {
    return nmax;
  }
  
  //! Get number of chain members
  /*! \return Number of chain members in particle address list #p
   */
  std::size_t getNchain() const {
    return nchain;
  }

  //! Get particle masses and store into array
  /*! Obtain particle masses and store into the double array
    @param[in] mass: double array to store the particle masses
   */
  void getMassAll(double mass[]) {
    for (std::size_t i=0;i<num;i++) mass[i] = (*this)[i].getMass();
  }

  //! Add new particle
  /*! Add new particle address at the end of the particle address list #p
    @param [in] a: new particle
   */
  void add(particle &a) {
    if (alloc_flag) {
      std::cerr<<"Error: chainlist already allocate memory for particle list, to avoid confusion, no new particle address can be appended!\n";
      abort();
    }
    if (num<nmax) {
      cflag[num] = false;
      p[num] = &a;
      num++;
    }
    else {
      std::cerr<<"Error: chainlist overflow! maximum number is "<<nmax<<std::endl;
      abort();
    }
  }

  //! Add a list of particles
  /*! ADD a list of particles at the end of the particle address list #p
    @param [in] n: number of new particles
    @param [in] a: array of new particles
   */
  void add(const std::size_t n, particle a[]) {
    if (alloc_flag) {
      std::cerr<<"Error: chainlist already allocate memory for particle list, to avoid confusion, no new particle address can be appended!\n";
      abort();
    }
    if (num+n<=nmax) {
      for (std::size_t i=0;i<n;i++) {
        cflag[num+i] = false;
        p[num+i] = &a[i];
      }
      num +=n;
    }
    else {
      std::cerr<<"Error: chainlist overflow! maximum number is "<<nmax<<", current number "<<num<<", try to add "<<n<<std::endl;
      abort();
    }
  }

  //! Add a chain in particle address list #p
  /*! Add one chain's address at the end of particle address list #p
    @param [in] a: new chain particle
   */
  void add(chain<particle,int_par> &a) {
    if (alloc_flag) {
      std::cerr<<"Error: chainlist already allocate memory for particle list, to avoid confusion, no new particle address can be appended!\n";
      abort();
    }
    if (num<nmax) {
      cflag[num] = true;
      p[num] = &a;
      num++;
      nchain++;
    }
    else {
      std::cerr<<"Error: chainlist overflow! maximum number is "<<nmax<<std::endl;
      abort();
    }
  }

  //! remove one particle
  /*! remove one particle from #p
     @param [in] i: particle index in #p needs to be removed
     @param [in] option: position update option: 
                 - true: shift last particle to current position (defaulted);
                 - false: shift all right particle to left by one
  */
  void remove(const std::size_t i, bool option=true) {
    if (option) {
      if (i<num-1) {
        if (cflag[i]) nchain--;
        num--;
        cflag[i] = cflag[num];
        p[i] = p[num];
        if (alloc_flag) p[num]=NULL;
      }
      else {
        std::cerr<<"Warning!: try to remove non-existing particle (index "<<i<<" > maximum number "<<num<<"std::endl";
      }
    }
    else {
      if (i<num-1) {
        if (cflag[i]) nchain--;
        num--;
        for (std::size_t j=i;j<num;j++) {
          cflag[j] = cflag[j+1];
          p[j] = p[j+1];
        }
        if (alloc_flag) p[num]=NULL;
      }
      else {
        std::cerr<<"Warning!: try to remove non-existing particle (index "<<i<<" > maximum number "<<num<<"std::endl";
      }
    }
  }

  //! Return the i^th of memebr's reference in the list
  /*! [] Operator overloading, return the i^th particle reference from the particle address list #p
    @param [in] i: the index of member in the list #p
    \return 
           - if \a i^th member is particle, return particle reference
           - if \a i^th member is a chain, return the center-of-mass particle in the chain (\ref ARC::chain.cm)
   */
  particle &operator [](const std::size_t i) const {
    if (i>=num) {
      std::cerr<<"Error: the required index "<<i<<" exceed the current particle list boundary (total number = "<<num<<std::endl;
      abort();
    }
    if (cflag[i]) {
      return ((chain<particle,int_par>*)p[i])->cm;
    }
    else {
      return *((particle*)p[i]);
    }
  }

  //! Dump function for particle data
  /*! Dump particle data into file using fwrite (binary format). Number of particles is written first, then the data of particle class in the particle list #p.
    @param [in] filename: file to store the data
   */
  void dump(const char* filename) {
    std::FILE* pout = std::fopen(filename,"w");
    if (pout==NULL) std::cerr<<"Error: filename "<<filename<<" cannot be open!\n";
    else {
      fwrite(&num,sizeof(int),1,pout);
      for (std::size_t i=0; i<num; i++) {
        fwrite(&((*this)[i]),sizeof(particle),1,pout);
      }
      fclose(pout);
    }
  }

  //! Read function for particle data
  /*! Read particle data using fread (binary format). Number of particles should be first variable, then the data of particles (the data structure should be consistent with particle class).
      Notice if read is used, the particle data memory is allocated inside chainlist, this is different from add() function.
    @param [in] filename: file to read the data
   */
  void read(const char* filename) {
    std::FILE* pin = std::fopen(filename,"r");
    if (pin==NULL) std::cerr<<"Error: filename "<<filename<<" cannot be open!\n";
    else {
      int n;
      std::size_t rn = fread(&n,sizeof(int),1,pin);
      if(rn<1) {
        std::cerr<<"Error: cannot read particle number!\n";
        abort();
      }

      allocate(n);
      for (std::size_t i=0; i<n; i++) {
        rn = fread(p[i+num],sizeof(particle),1,pin);
        if (rn<1) {
          std::cerr<<"Error: cannot read particle "<<i<<", total particle number is "<<n<<"!\n";
          abort();
        }
      }
      fclose(pin);
    }
  }

  //! Is i^th a chain?
  /*! Check whether the i^th member in the particle address list is chain
     @param [in] i: member index in list #p
     \return  true: chain; false: particle
   */
  bool isChain (const std::size_t i) const {
    return cflag[i];
  }
  
  //! Get particle list in a chain type member
  /*! Return the address of particle list of the i^th member in particle #p (\ref ARC::chain.p)
    @param [in] i: member index in list the particle address list #p
    \return
     - if i^th member in #p is chain, returen the address of particle list (\ref ARC::chain.p).
     - else return NULL
   */
  chain<particle, int_par> *getSub (const std::size_t i) const {
    if (cflag[i]) return (chain<particle, int_par>*)p[i];
    else return NULL;
  }

};

}
