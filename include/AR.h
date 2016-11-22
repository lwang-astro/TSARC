#pragma once

#include <iostream>
#include <cstring>
#include <cmath>
#include <cstdlib>
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
struct timeprofile{
public:
  double t_apw;    ///< APW calculation
  double t_uplink; ///< update link
  double t_lf;     ///< leap-frog
  double t_ep;     ///< extrapolation
  double t_pext;   ///< perturber force

  int* stepcount;  ///< iteration step count in extrapolation

  timeprofile() {reset_tp();}  ///< initialization 

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
    stepcount=NULL;
  }

  ~timeprofile() { if (stepcount) delete[] stepcount;}

};

#endif

//declaration
template <class particle> class chain;
template <class particle> class chainlist;
class chainpars;

//! Function pointer type to function for calculation of acceleration (potential) and time transformation function dW/dr (W) from particle j to particle i.
/*!          @param[out] A: acceleration vector.
             @param[out] P: potential (positive value).
             @param[out] dW: dW/dr (used for beta>0) .
             @param[out] W: time transformation function (used for beta>0).
             @param[in] X: relative position (1:3).
             @param[in] mi: particle i mass.
             @param[in] mj: particle j mass.
             @param[in] mm2: smooth mass coefficient \f$\sum m_i  m_j /(N (N-1)/2) (i<j) \f$ (Notice it is only calculated when m_smooth is true in chainpars).
             @param[in] epi: adjustable parameter (set in ARC::chainpars).
 */
typedef void (*pair_AW) (double *, double &, double *, double &, const double *, const double &, const double &, const double &, const double &);

//! Function pointer type to function for calculation of acceleration (potential) from particle j to particle i.
/*!         @param[out]  A: acceleration vector.
            @param[out]  P: potential (positive value).
            @param[in]  pi: position vector i.
            @param[in]  pj: position vector j.
            @param[in]  mi: particle mass i.
            @param[in]  mj: particle mass j.
 */
typedef void (*pair_Ap) (double *, double &, const double*, const double*, const double&, const double&);

//! Newtonian acceleration and dW/dr from particle k to particle j (function type of \link #ARC::pair_AW \endlink) 
/*!          @param[out] A: Newtonian acceleration vector. \f$A[1:3] = m_j m_k xjk[1:3] / |xjk|^3 \f$.
             @param[out] P: Newtonian potential. \f$ P = m_j m_k /|xjk| \f$
             @param[out] dW: TTL time transformation function derivates based on position vector dW/dr (used for TTL method). \f$dW[1:3] = w_{jk} xjk[1:3] /|xjk|^3 \f$
             @param[out] W: TTL time transformation function (used for TTL method) \f$ W = \sum_{j<k} w_{jk} /|xjk| \f$
             @param[in] xjk: relative position vector [1:3] from particle j to particle k
             @param[in] mj: particle j mass.
             @param[in] mk: particle k mass.
             @param[in] mm2: smooth mass coefficient \f$ \sum_{j<k} m_j m_k /(N(N-1)/2) \f$ (Notice in \ref ARC::chain, it is only calculated when ARC::chainpars::m_smooth is true).
             @param[in] epi: adjustable parameter (ARC::chainpars::m_epi).
                             If epi>0:  \f$ w_{jk} = mm2 \f$ (\f$ mm2 = \sum_{j<k} m_j m_k / (N( N - 1 )/2) \f$) else: \f$w_{jk} = m_j m_k\f$.\n
*/
void Newtonian_AW (double A[3], double &P, double dW[3], double &W, const double xjk[3], const double &mj, const double &mk, const double &mm2, const double &epi) {

  // distance
  double rjk = std::sqrt(xjk[0]*xjk[0]+xjk[1]*xjk[1]+xjk[2]*xjk[2]);  

  // mass parameters
  double mmjk = mj*mk; // m_i*m_j
  double wjk;
  if (mm2>0 && epi>0) {
    // Wjk = mm2 if m_i*m_j < epi*m'^2; 0 otherwise;
    if (mmjk<epi*mm2) wjk = mm2;
    else wjk = 0;
  }
  else {
    wjk = mmjk;    // Wjk = m_i*m_j
  }
  
  P = mmjk / rjk;  // Potential energy
  W = wjk / rjk;   // Transformation coefficient
        
  // Acceleration
  double rjk3 = rjk*rjk*rjk;
  double mor3 = mk / rjk3;
  A[0] = mor3 * xjk[0];
  A[1] = mor3 * xjk[1];
  A[2] = mor3 * xjk[2];

  // dW/dr
  mor3 = wjk / rjk3;
  dW[0] = mor3 * xjk[0];
  dW[1] = mor3 * xjk[1];
  dW[2] = mor3 * xjk[2];
  
}

//! Newtonian acceleration from particle p to particle i (function type of ::ARC::pair_Ap)
/*! 
  @param[out]  A: acceleration vector. \f$A[1:3] = m_i m_p (xp[1:3]-xi[1:3]) / |xp-xi|^3 \f$.
  @param[out]  P: potential. \f$ P = m_i m_p /|xp-xi|^3\f$
  @param[in]  xi: position vector i.
  @param[in]  xp: position vector p.
  @param[in]  mi: particle mass i.
  @param[in]  mp: particle mass p.
 */
void Newtonian_Ap (double A[3], double &P, const double xi[3], const double xp[3], const double &mi, const double &mp){
  double dx = xp[0] - xi[0];
  double dy = xp[1] - xi[1];
  double dz = xp[2] - xi[2];

  double dr2 = dx*dx + dy*dy + dz*dz;
  double dr  = std::sqrt(dr2);
  double dr3 = dr*dr2;

  A[0] = mp * dx / dr3;
  A[1] = mp * dy / dr3;
  A[2] = mp * dz / dr3;

  P = mi*mp / dr;
  
}
 
//! The chain parameter controller class
/*!
  This class control the acceleration function (::ARC::pair_AW, ::ARC::pair_Ap), integration methods and error parameters for \ref chain.
 */
class chainpars{
template <class T> friend class chain;
private:
  // pair force
  pair_AW pp_AW;  ///< acceleration and dW/dr of two particle
  pair_Ap pp_Ap;  ///< accelaration of two particle
  
  // time step integration parameter
  double alpha; ///< logH cofficient
  double beta;  ///< TTL cofficient
  double gamma; ///< constant

  // mass coefficients parameter
  double m_epi; ///< smooth parameter
  bool m_smooth; ///< whether to use smooth mass coefficients

  // time step
  double dtmin; ///< minimum physical time step
  double dterr; ///< physical time error criterion

  // extrapolation control parameter
  double exp_error;        ///< relative error requirement for extrapolation
  std::size_t exp_itermax; ///< maximum times for iteration.
  int exp_method;         ///< 1: Romberg method; others: Rational interpolation method
  int exp_sequence;       ///< 1: even sequence {h, h/2, h/4, h/8 ...}; 2: Bulirsch & Stoer sequence {h, h/2, h/3, h/4, h/6, h/8 ...}; other. 4k sequence {h/2, h/6, h/10, h/14 ...}

  int* step; ///< substep sequence
  std::size_t  opt_iter; ///< optimized iteration index

  int** bin_index; ///< binomial coefficients

public:

  //! constructor with defaulted parameters
  /*! - Acceleration function use ARC::Newtonian_AW() and ARC::Newtonian_Ap().
      - ARC method use logarithmic Hamiltonian (logH) (#alpha = 1.0, #beta = 0.0 #gamma = 0.0).
      - Smooth mass cofficients is used (#m_smooth = true) with smooth parameter #m_epi = 0.001.
      - Phase/energy error limit #exp_error = 1e-10.
      - Minimum physical time step #dtmin = 5.4e-20.
      - Time synchronization error limit #dterr = 1e-6.
      - Maximum extrapolation iteration number #exp_itermax = 20
      - Bulirsch & Stoer sequence {h, h/2, h/3, h/4, h/6...} is used
      - Optimized extrapolation interation order for auto-adjust integration step size #opt_iter = 5
   */
  chainpars(): alpha(1.0), beta(0.0), gamma(0.0), m_epi(0.001), m_smooth(true) {
    step = NULL;
    bin_index = NULL;
    setEXP(1E-10, 5.4E-20, 1E-6, 20, 2, 2, 5);
    pp_AW = &Newtonian_AW;
    pp_Ap = &Newtonian_Ap;
  }

  //! constructor
  /*!
    @param [in] aw: acceleration and dW/dr of two particle calculation function pointer with ::ARC::pair_AW type. (interaction between memebrs)
    @param [in] ap: acceleration calculation function pointer with ::ARC::pair_Ap type. (interaction between member and perturber)
    @param [in] a,b,g: ARC time transformation method coefficients (\f$ dt = ds/[a *(logH) + b * (TTL) + g])\f$. \n
                - a: Logarithmic Hamiltonian (logH) method coefficient (0.0, 1.0)
                - b: Time-Transformed Leapfrog (TTL) method coefficient (0.0, 1.0)
                - g: Constant coefficient (no time transformation) 
    @param [in] eps: Smooth mass coefficient parameter (defaulted #m_epi = 0.001)
    @param [in] mm: Whether to use smooth mass coefficients (defaulted #m_smooth = true)
    @param [in] error: Phase/energy error limit (defaulted #exp_error = 1e-10)
    @param [in] dtm: Minimum physical time step (defaulted #dtmin = 5.4e-20)
    @param [in] dte: Time synchronization error limit (defaulted #dterr = 1e-6)
    @param [in] itermax: Maximum extrapolation iteration number (defaulted #exp_itermax = 20)
    @param [in] ext_method: 1: Romberg interpolation method; others: Rational interpolation method (defaulted: Rational)
    @param [in] ext_sequence: 1: even sequence {h, h/2, h/4, h/8 ...}; 2: Bulirsch & Stoer (BS) sequence {h, h/2, h/3, h/4, h/6, h/8 ...}; 3: 4k sequence {h, h/2, h/6, h/10, h/14 ...}; others: Harmonic sequence {h, h/2, h/3, h/4 ...} (defaulted 2. BS sequence)
    @param [in] optiter: Optimized extrapolation interation order for auto-adjust integration step size (defaulted #opt_iter = 5)
   */
  chainpars(pair_AW aw, pair_Ap ap, const double a, const double b, const double g, const double eps=0.001, const bool mm=true, const double error=1E-10, const double dtm=5.4e-20, const double dte=1e-6, const std::size_t itermax=20, const int ext_method=2, const int ext_sequence=2, const std::size_t optiter=5) {
    step = NULL;
    bin_index = NULL;
    setabg(a,b,g);
    setM(eps,mm);
    setEXP(error,dtm,dte,itermax,ext_method,ext_sequence,optiter);
    setA(aw,ap);
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

  //! Set acceleration (dW/dr) function
  /*! Set pair acceleration (dW/dr) function for interaction between chain members and from perturbers
    @param [in] aw: acceleration and dW/dr of two particle calculation function pointer with ::ARC::pair_AW type. (interaction between memebrs)
    @param [in] ap: acceleration calculation function pointer with ::ARC::pair_Ap type. (interaction between member and perturber)
  */
  void setA(pair_AW aw, pair_Ap ap) {
    pp_AW = aw;
    pp_Ap = ap;
    // safety check
    if (pp_AW==NULL) {
      std::cerr<<"Error: accelaration / W calculator function is NULL\n";
      abort();
    }
    if (pp_Ap==NULL) {
      std::cerr<<"Error: perturber accelaration calculator function is NULL\n";
      abort();
    }
  }
  
     
  //! Set time step transformation parameter
  /*!
      Set parameter a,b,g where physical time step \f$ dt = ds/[a *(logH) + b * (TTL) + g]\f$ \n
      @param [in] a: Logarithmic Hamiltonian (logH) method coefficient (0.0, 1.0)
      @param [in] b: Time-Transformed Leapfrog (TTL) method coefficient (0.0, 1.0)
      @param [in] g: Constant coefficient (no time transformation) 
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
    if (alpha==0&&gamma==0&&m_smooth) {
      std::cerr<<"TTL method is used with gamma=0, thus smooth mass coefficient is forced to be switched off\n";
      setM(0.001,false);
    }
  }

  //! Set smooth mass coefficient parameter
  /*! Set smooth mass coefficient parameter epi and determine whether use this (#m_smooth) \n
      The smooth mass coefficient is defined as \f$ mm2 = \sum m_i  m_j /(N (N-1)/2) (i<j) \f$.
     @param[in] epi: smooth parameter 
     @param[in] mm: m_smooth (whether use smooth mass coefficients)
  */
  void setM(const double epi=0.001, const bool mm=true) {
    m_epi = epi;
    m_smooth = mm;
    if (m_epi==0&&m_smooth) {
      std::cerr<<"Error: smooth mass coefficients are used, but smooth coefficient epi is zero!\n";
      abort();
    }
    if (alpha==0&&gamma==0) {
      std::cerr<<"Warning: alpha=0 and gamma=0, the smooth mass coefficients may cause initial w = 0 and result in zero time step!";
    }
  }

  //! Set extrapolation parameters
  /*! Set extrapolation related parameters as following:
    @param [in] error: phase/energy relative error requirement for extrapolation (defaulted #exp_error = 1e-10)
    @param [in] dtm: minimum physical time step allown (defaulted #dtmin = 5.4e-20)
    @param [in] dte: Time synchronization error limit (defaulted #dterr = 1e-6)
    @param [in] itermax: maximum order (index in sequence) for extrapolation iteration. (defaulted #exp_itermax = 20)
    @param [in] methods: 1: Romberg method; others: Rational interpolation method (defaulted Rational)
    @param [in] sequences: 1: even sequence {h, h/2, h/4, h/8 ...}; 2: Bulirsch & Stoer (BS) sequence {h, h/2, h/3, h/4, h/6, h/8 ...}; 3: 4k sequence {h, h/2, h/6, h/10, h/14 ...}; others: Harmonic sequence {h, h/2, h/3, h/4 ...} (defaulted 2. BS sequence)
    @param [in] optiter: Optimized interation order for auto-adjust integration time step (defaulted #opt_iter = 5)
  */
  void setEXP(const double error=1E-10, const double dtm=5.4e-20, const double dte=1e-6, const std::size_t itermax=20, const int methods=2, const int sequences=2, const std::size_t optiter=5) {
    exp_error = error;
    exp_method = methods;
    exp_sequence = sequences;
    dterr = dte;
    dtmin = dtm;
    opt_iter=optiter;

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
    if (sequences==1) EP::seq_Romberg(step,itermax+1);
    // Bulirsch & Stoer sequence {h, h/2, h/3, h/4, h/6, h/8 ...}
    else if (sequences==2) EP::seq_BS(step,itermax+1);
    // E. Hairer (4k) sequences {h, h/2, h/6, h/10, h/14 ...}
    else if (sequences==3) EP::seq_Hairer(step,itermax+1);
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
  
};

//! ARC class based on template class particle
/*!
  Major class for ARC integration of few bodies
  
  It depend on the template class particle. This particle class should contain public member functions for reading and writing mass, position and velocity (see sample in Particle::setPos(), Particle::setVel(), Particle::setMass(), Particle::getPos(), Particle::getVel(), Particle::getMass())

  The basic way to use ARC integration is shown as following:
  1. Construct a chain class with template class particle and a parameter controller of \ref ARC::chainpars. (The \ref ARC::chainpars should be configured first before doing integration. see its document for detail).
  2. Add existed particle 'A' (or a list of particles, or a chain type particle) into chain particle list (\ref chain.p) using chain.addP(). Notice the chain.addP() only registers the particle A's memory address into \ref chain.p without copying data. The chain integration will directly modify the position and velocity of particle A.
  3. Add perturbers into chain perturber list (\ref chain.pext) using chain.addPext() (also only register the particle address)
  4. Initialization chain with chain.init(). Notice this function is necessary to be called before integration. Also be careful that after this initialization, the positions and velocites of particles registered in \ref chain.p will be shifted from their original frame to their center-of-mass frame. The particle type member variable \ref chain.cm stores the center-of-mass data of these particles (the mass of \ref chain.cm is the total mass of all member particles).
  5. Call integration functions (chain.Leapfrog_step_forward() or chain.extrapolation_integration()). The former use only Leapfrog method and the latter use extrapolation method to obtain high accuracy of integration.
  6. After call integration functions, the particles are integrated to new time. Because in ARC method, the time is also integrated and cannot be predicted before integration, thus the iteration need to be done to get correct physical time you want (see detailed in chain.extrapolation_integration() document).
  7. Notice that after chain.init() and integration, the particles are always in the center-of-mass frame. If you want to shift them back to the original frame, the chain.center_shift_inverse() should be used. But after this function is used, you should use chain.center_shift() before the next integration.
 */
template <class particle>
class chain{
  typedef double double3[3];
  double3 *X;  ///< relative position
  double3 *V;  ///< relative velocity
  std::size_t *list;   ///< chain index list
  double3 *acc; ///< acceleration
  double3 *pf;  ///< perturber force
  double3 *dWdr; ///< \partial Omega/ \partial rk

  //integration parameters=======================================//
  double t;    ///< time
  double w;    ///< time transformation parameter
  double B;    ///< Binding energy (time momentum)

  //template parameters==========================================//
  double W;     ///< time transformation function
  double Ekin;  ///< kinetic energy
  double Pot;   ///< potential
  double mm2;   ///< mean mass production \sum m_i*m_j/(N(N-1)/2) (i<j)

  //number =======================================================//
  std::size_t num;      ///< total number of chain particles
  std::size_t nmax;     ///< maximum number 

  //monitor flags
  bool F_Pmod;     ///< indicate whether particle list is modified (true: modified)
  int  F_Porigin;  ///< indicate whether particle is shifted back to original frame (1: original frame: 0: center-of-mass frame; 2: only position is original frame)

  const chainpars *pars;   ///< chain parameter controller

public:

  particle cm;              ///< center mass particle
  chainlist<particle> p;    ///< particle list
  chainlist<particle> pext; ///< perturber list

#ifdef TIME_PROFILE
  timeprofile profile;
#endif
  
  //! Constructor
  /*! Construct chain with allocated memory
      @param [in] n: maximum number of particles (will be used to allocate memory)
      @param [in] par: chain option controller class \ref ARC::chainpars
   */
  chain(std::size_t n, const chainpars &par):  pars(&par) {
    nmax=0;
    mm2=-std::numeric_limits<double>::infinity();
    allocate(n);
  }

  //! Constructor
  /*! Construct chain without memory allocate, need to call allocate() later. 
     @param [in] par: chain option controller class \ref ARC::chainpars
   */
  chain(const chainpars &par): pars(&par), F_Pmod(false), F_Porigin(1), num(0), nmax(0) {
    mm2=-std::numeric_limits<double>::infinity();
  }

  //! Allocate memory
  /*! Allocate memory for maximum particle number n
     @param [in] n: maximum number of particles
   */
  void allocate(std::size_t n) {
    if (nmax) {
      std::cerr<<"Error: chain memory allocation is already done\n";
      abort();
    }
    num = n;
    nmax = n;
    X=new double3[n-1];
    V=new double3[n-1];
    list=new std::size_t[n];
    acc=new double3[n];
    pf=new double3[n];
    dWdr=new double3[n];
    p.init(n);
    F_Pmod=false;
    F_Porigin=1;
    mm2=-std::numeric_limits<double>::infinity();
  }

  //! Clear function
  /*! Clear allocated memory and set maximum number of particle to zero
   */
  void clear() {
    if (nmax>0) {
      delete[] X;
      delete[] V;
      delete[] list;
      delete[] acc;
      delete[] pf;
      delete[] dWdr;
      num = 0;
      nmax = 0;
    }
    F_Pmod=false;
    F_Porigin=1;
    mm2=-std::numeric_limits<double>::infinity();
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
  }

private:
  //! Update number of particles
  /*! Update the number of particle (#num) due to current particle list number
     @param [in] n: current particle number in particle list #p
  */
  void update_num(const std::size_t n) {
    if (n>nmax) {
      std::cerr<<"Error: particle number "<<n<<" is larger than Chain number limit "<<num<<std::endl;
      abort();
    }
    else{
      num = n; ///< Update Chain number
    }
  }

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


  //! Calculate smooth particle coefficient
  /*! Get smooth particle coefficient \f$ mm2 = \sum m_i  m_j /(N (N-1)/2) (i<j) \f$ based on masses of particles in #p
   */
  void calc_mm2() {
    // calcualte m'^2
    mm2 = 0;
    for (std::size_t i=0;i<num;i++) {
      for (std::size_t j=i+1;j<num;j++) {
        mm2 += p[i].getMass() * p[j].getMass();
      }
    }
    mm2 /= num * (num - 1) / 2;
  }    

  
  //! Calculate acceleration (potential) and transformation parameter
  /*! Get distance Matrix, acceleration, dm/dr, potential and transformation parameter
      based on particle masses in #p, using current #X, #V
      (notice the acceleration and dW/dr array index follow particle p to avoid additional shift when chain list change).

      @param [in] force: external force (acceleration) for each particle, (not perturber forces)
      @param [in] resolve_flag: flag to determine whether to resolve sub-chain particles for force calculations. (defaulted false)
  */
  void calc_rAPW (const double3 *force=NULL, const bool resolve_flag=false) {
#ifdef TIME_PROFILE
    profile.t_apw -= get_wtime();
#endif
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
        pars->pp_AW(At, Pt, dWt, Wt, xjk, mj, mk, mm2, pars->m_epi);

        // resolve sub-chain
        if(resolve_flag && p.isChain(lk)) {
          chain<particle>*ck = p.getSub(lk);
          // center shift to current frame
          ck->center_shift_inverse_X();
          const std::size_t cn = ck->p.getN();
          Pt = 0;
          for (std::size_t i=0;i<3;i++) At[i]=0.0;
          for (std::size_t i=0;i<cn;i++) {
            double Ptemp;
            double3 Atemp;
            pars->pp_Ap(Atemp, Ptemp, xj, ck->p[i].getPos(), mj, ck->p[i].getMass());

            // Acceleration
            At[0] += Atemp[0];
            At[1] += Atemp[1];
            At[2] += Atemp[2];
            
            // Potential
            if (k>j) Pt += Ptemp;
          }
          // center shift back
          ck->center_shift_X();
        }

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
      // add external acceleration
      if (force!=NULL) {
        acc[lj][0] += pf[lj][0] + force[lj][0];
        acc[lj][1] += pf[lj][1] + force[lj][1];
        acc[lj][2] += pf[lj][2] + force[lj][2];
      }
      else {
        acc[lj][0] += pf[lj][0];
        acc[lj][1] += pf[lj][1];
        acc[lj][2] += pf[lj][2];
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

  //! Calculate physical time step for X
  /*! Calculate physical time step dt for #X based on ds
     @param [in] ds:  integration step size (not physical time step) 
     \return     dt:  physical integration time step for #X
  */
  double calc_dt_X(const double ds) {
    // determine the physical time step
    double dt = ds / (pars->alpha * (Ekin + B) + pars->beta * w + pars->gamma);
    if (std::abs(dt) < pars->dtmin) {
      std::cerr<<"Warning!: physical time step too small: "<<dt<<std::endl;
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
    return ds / (pars->alpha * Pot + pars->beta * W + pars->gamma);
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

  //! Step forward of B and w
  /*! One step integration of B and w.
      - \f$B += dt * \sum ( - m_k * <v_k> \dot f_k)\f$ 
      - \f$w += dt * \sum ( dm/dr_k \dot v_k)\f$ 
     @param [in] dt: time step for V
     @param [in] ave_v: averaged velocity
     @param [in] force: external force
     @param [in] p: particle list (only use mass)
     @param [in] fpf: if perturber force is not zero, it is true
   */
  void step_forward_Bw(const double dt, const double3* ave_v, const double3* force, const bool fpf) {
    double dB = 0.0;
    double dw = 0.0;
    if (force!=NULL||fpf||pars->beta>0) {
      for (std::size_t i=0;i<num;i++) {
        if (force!=NULL) {
          dB -= p[i].getMass() * ( ave_v[i][0] * (pf[i][0] + force[i][0]) 
                            + ave_v[i][1] * (pf[i][1] + force[i][1]) 
                            + ave_v[i][2] * (pf[i][2] + force[i][2]));
        }
        else if (fpf){
          dB -= p[i].getMass() * ( ave_v[i][0] * pf[i][0] 
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
    B += dt * dB;
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
    // resolve current V
    double3 vc={0};
    double3 xc={0};
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
      if (i==0) {
        const double mk = p[lk].getMass();
        xc[0] += mk * rk[0];
        xc[1] += mk * rk[1];
        xc[2] += mk * rk[2];
        vc[0] += mk * vk[0];
        vc[1] += mk * vk[1];
        vc[2] += mk * vk[2];
      }        
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
    double3 xc={0};
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
      if (i==0) {
        const double mk = p[lk].getMass();
        xc[0] += mk * rk[0];
        xc[1] += mk * rk[1];
        xc[2] += mk * rk[2];
      }
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
    double3 vc={0};
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
      if (i==0) {
        const double mk = p[lk].getMass();
        vc[0] += mk * vk[0];
        vc[1] += mk * vk[1];
        vc[2] += mk * vk[2];
      }        
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
            chain<particle>*cj = pext.getSub(j);
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
              pars->pp_Ap(Atemp, Pt, xi, xk, mi, cj->p[k].getMass());

              // Acceleration
              At[0] += Atemp[0];
              At[1] += Atemp[1];
              At[2] += Atemp[2];
            }
          }
          else {
            // perturber force
            pars->pp_Ap(At, Pt, xi, pext[j].getPos(), mi, pext[j].getMass());
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
    std::memcpy(Xbk,X,(num-1)*3*sizeof(double));
    
    // backup current V
    double3* Vbk = new double3[num-1];
    std::memcpy(Vbk,V,(num-1)*3*sizeof(double));

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
    double cmr[3]={};
    double cmv[3]={};
    double cmm = 0;
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
      for (std::size_t i=0;i<num;i++) {
        const double *ri = p[i].getPos();
        const double *rc = cm.getPos();
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
      for (std::size_t i=0;i<num;i++) {
        const double *ri = p[i].getPos();
        const double *rc = cm.getPos();
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
        double dts=1/(pars->alpha * (Ekin + B) + pars->beta * w + pars->gamma);
        
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
#ifdef DEBUG
              std::cerr<<"Poly_coff: n= "<<n<<"; i="<<i<<"; ik="<<ik<<"; coff="<<coff<<std::endl;
#endif
              dpoly[j][0] += coff * dts;
//              if (!tflag) {
//                dpoly[j][1] += coff * B;
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
        double dts=1.0/(pars->alpha * (Ekin + B) + pars->beta * w + pars->gamma);

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
//              dpoly[j][1] += coff * B;
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
//              dpoly[j][1+ir] += coff * B;
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
  void addP(chain<particle> &a) {
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


  //! Allocate memory for perturber list
  /*! Allocate memory for perturber particle list with maximum number of \a n
     @param [in] n: maximum number of perturbers
   */
  void initPext(const std::size_t n) {
    if (pext.getN()) {
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
  void addPext(chain<particle> &a) { pext.add(a);}

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
  
  //! Initialization
  /*! Initialize chain based on particle list #p. After this function, positions and velocities of particles in #p will be shifted to their center-of-mass frame. \n
      Chain order list #list, relative position #X, velocity #V, initial system energy #B and initial time transformation parameter #w are calculated.
      The particle modification indicator (isPmod()) will be set to false.
    @param [in] time: current time of particle system
    @param [in] force: external force
  */
  void init(const double time, const double3* force=NULL) {
    // update number indicator
    update_num(p.getN());
    
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
    
    // if smooth mass coefficients are used, calculate mm2;
    if (pars->m_smooth) calc_mm2();

    // set relative distance matrix, acceleration, potential and transformation parameter
    calc_rAPW(force);

    // Initial intgrt value t
    t = time;

    // kinetic energy
    calc_Ekin();

    // initial time step parameter
    B = Pot - Ekin;
    w = W;

    // set F_Pmod to false
    F_Pmod = false;

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
    std::cerr<<std::setw(WIDTH)<<Ekin<<std::setw(WIDTH)<<Pot<<std::setw(WIDTH)<<B+Ekin-Pot<<std::setw(WIDTH)<<B<<std::setw(WIDTH)<<w<<std::endl;
    #endif*/
    
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
      for (std::size_t i=0;i<num;i++) {
        if (F_Porigin==0) {
          const double *ri = p[i].getPos();
          const double *rc = cm.getPos();
          p[i].setPos(ri[0] + rc[0],
                      ri[1] + rc[1],
                      ri[2] + rc[2]);
        }
        const double *vi = p[i].getVel();
        const double *vc = cm.getVel();
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
      for (std::size_t i=0;i<num;i++) {
        const double *ri = p[i].getPos();
        const double *rc = cm.getPos();
        p[i].setPos(ri[0] - rc[0],
                    ri[1] - rc[1],
                    ri[2] - rc[2]);
        if (F_Porigin==1) {
          const double *vi = p[i].getVel();
          const double *vc = cm.getVel();
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

  //! Backup chain data (#t, #B, #w, #X, #V)
  /*! Backup chain data to one dimensional array. 
      - #t : current physical time
      - #B : current system energy
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
    db[1] = B;
    db[2] = w;
    std::memcpy(&db[3], X, ndata*sizeof(double));
    std::memcpy(&db[3+ndata], V, ndata*sizeof(double));
  }
     
  //! Restore chain data (#t, #B, #w, #X, #V)
  /*! Restore integration data from one dimensional array, the order of data should be #t, #B, #w, #X[#num][3], #V[#num][3]
      - #t : current physical time
      - #B : current system energy
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
    B = db[1];
    w = db[2];
    std::memcpy(X, &db[3], ndata*sizeof(double));
    std::memcpy(V, &db[3+ndata], ndata*sizeof(double));
  }

  //! Leapfrog integrator
  /*! Integration with Leapfrog method. \n
      The positions and velocities of particles in #p will be integrated in the center-of-mass frame
      @param [in] s: Integration step size
      @param [in] n: number of sub-steps needed to be divided. Integration step size is (\a s/\a n) and do \a n times as: X(s/2n)V(s/n)X(s/n)V(s/n)..X(s/2n)
      @param [in] force: external force (not perturber forces which are calculated in pert_force)
      @param [in] check_flag: 2: check link every step; 1: check link at then end; 0 :no check
      @param [in] dpoly: two dimensional array for storing 0 to \a (ndmax-1)'th order central difference of physical time #t at \a s/2 as a function of \a s, array size should be [ndmax][1]
      @param [in] ndmax: dpoly array size and the maximum difference is ndmax-1
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
    if (pars->m_smooth&&mm2==-std::numeric_limits<double>::infinity()) {
      std::cerr<<"Error: smooth mass coefficients are used, but averaged mass coefficient mm2 is not calculated!\n";
      abort();
    }

    double ds = s/double(n);
    double3* ave_v=new double3[num];  // average velocity
    bool fpf = false;                 // perturber force indicator
    const int np = pext.getN();
    if (np>0) fpf = true;

    // for polynomial coefficient calculation first point
    // middle difference (first array is used to store the f_1/2)
    if(pars->exp_sequence==3) mid_diff_calc(&dpoly[2],ndmax-2,0,n);
    // edge difference
    else {
      dpoly[0][0]=1.0/(pars->alpha * Pot + pars->beta * W + pars->gamma);
      edge_diff_calc(&dpoly[1],ndmax-1,0,n);
    }
                                               
    
    // integration loop
    for (std::size_t i=0;i<n;i++) {
      // half step forward for t (dependence: Ekin, B, w)
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
        chain<particle> **clist = new chain<particle>*[nct];

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

      // forward B and w (dependence: dt(V), ave_v, p.m, p.v, force, dWdr, pf)
      step_forward_Bw(dvt,ave_v,force,fpf);

      // Calcuale Kinetic energy (dependence: p.m, p.v)
      calc_Ekin();

      // step forward for t (dependence: Ekin, B, w)
      dt = calc_dt_X(ds*0.5);
      t += dt;

      // step forward for X (dependence: X, V)
      step_forward_X(dt);
      
      // for interpolation polynomial coefficient (difference)
      // middle difference (first array is used to store the f_1/2)
      if(pars->exp_sequence==3) {
        if (i==n/2-1) {
          dpoly[0][0]=t; // y
          dpoly[1][0]=2.0*dt/ds; // f(x)
#ifdef DEBUG
          std::cerr<<"Mid time = "<<t<<", n="<<n<<"; i="<<i+1<<std::endl;
#endif
        }
        mid_diff_calc(&dpoly[2],ndmax-2,i+1,n);
      }
      // edge difference
      else edge_diff_calc(&dpoly[1],ndmax-1,i+1,n);
    }

    // resolve X at last, update p.x (dependence: X)
    resolve_X();

    // Update rjk, A, Pot, dWdr, W (notice A will be incorrect since pf is not updated)
    calc_rAPW(force);

    if(pars->exp_sequence!=3) dpoly[0][1]=1.0/(pars->alpha * Pot + pars->beta * W + pars->gamma);

#ifdef DEBUG
    std::cerr<<"Ending time = "<<t<<", n="<<n<<std::endl;
#endif
//#ifdef DEBUG
//    std::cerr<<std::setw(WIDTH)<<t;
//    //      for (int i=0;i<num;i++) 
//    //        for (int k=0;k<3;k++) std::cerr<<p[i].pos[k]<<" ";
//    for (int i=0;i<num-1;i++) {
//      for (int k=0;k<3;k++) std::cerr<<std::setw(WIDTH)<<X[i][k];
//      for (int k=0;k<3;k++) std::cerr<<std::setw(WIDTH)<<V[i][k];
//    }
//    std::cerr<<std::setw(WIDTH)<<Ekin<<std::setw(WIDTH)<<Pot<<std::setw(WIDTH)<<B+Ekin-Pot<<std::setw(WIDTH)<<B<<std::setw(WIDTH)<<w<<std::endl;
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
                      - if value is positive and after step \a ds, the ending physical time is larger than \a toff, the interpolation of physical time #t will be done instead of integration. In this case, positions and velocites of particles in #p and chain data are not changed, instead, returning value is the ds modification factor (negative value), which can be used to modified current \a ds and redo the integration by calling this function again with new ds to approach ending physical time of \a toff.
     @param [in] force: external force (not perturber forces which are calculated in pert_force)
     \return factor
            - if factor is positive, it is optimized step size modification factor for next step (\a ds *= factor)
            - if factor is negative and \a toff>0; the -factor is used for calculate new ds' = -factor * \a ds. Thus this function should be called again with new step size ds'. Thus the new result should have ending physical time close to \a toff.
   */
  double extrapolation_integration(const double ds, const double toff=-1.0, const double3* force=NULL) {
    // get parameters
    const double error = pars->exp_error;
    const std::size_t itermax = pars->exp_itermax;
    const int method = pars->exp_method;
    const int *step = pars->step;
    
#ifdef TIME_PROFILE
    profile.t_ep -= get_wtime();
#endif
    // array size indicator for relative position and velocity
    const std::size_t nrel = num-1;
    // data storage size (extra one is used to show the size of the array for safety check)
    const std::size_t dsize = 6*nrel+3;
    // edge difference array size;
    // const std::size_t psize=dsize*2;
    
    // for storage
    // for convenient, the data are storaged in one array with size (2*nrel+1)*3, which represents t, B, w, X[3][nrel], V[3][nrel]
    double d0[dsize+1],dtemp[dsize+1];
    double* dn[itermax];
    d0[dsize] = (double)dsize;    // label for safety check
    dtemp[dsize] = (double)dsize; // label for safety check
    for (std::size_t i=0; i<itermax; i++) {
      dn[i] = new double[dsize+1];
      dn[i][dsize] = (double)dsize; // label for safety check
    }
    double Ekin0;

    // for dense output polynomial
    bool ip_flag = true;  // interpolation coefficient calculation flag
    double** pd[itermax]; // central difference, [*] indicate different accuracy level 
    int ndmax[itermax];   // maximum difference order
    std::size_t pnn;      // data size
    if(pars->exp_sequence==3) pnn = 1;     // middle difference case
    else pnn = 2;         // edge two points case

    // for error check
    double3 CX,CXN;
    double cxerr=error+1.0;
    double eerr=error+1.0;
    double cxerr0=cxerr+1.0;
    double eerr0=eerr+1.0;
    // error for step estimation
    double werrmax=std::numeric_limits<double>::max();

    // backup initial values (t, B, w, X, V, Ekin)
    backup(d0);
    Ekin0 = Ekin;

    // new step
    double dsn = 1.0;
    
    std::size_t intcount = 0;  // iteration counter
    std::size_t itercount = 0; // iteration efforts count
    
    // first step
    
    while (cxerr > 0.5*error || std::abs(eerr) > 0.1*error) {
      // convergency check
      if (cxerr>=cxerr0 || std::abs(eerr) >= std::abs(eerr0)) {
        if (cxerr < error && std::abs(eerr) < error) break;
        else if (intcount > std::min((size_t)10,itermax)){
          std::cerr<<"Warning: extrapolation cannot converge anymore, energy error - current: "<<eerr<<"  previous: "<<eerr0<<"   , phase error - current: "<<cxerr<<"  previous: "<<cxerr0<<", try to change the error criterion (notice energy error is cumulative value)\n";
          break;
        }
      }
      // if completely converged, check energy error
      if (cxerr==0) {
        if (std::abs(eerr) > error)
          std::cerr<<"Warning: phase error reach zero but energy error "<<eerr<<" cannot reach criterion "<<error<<"!\n";
        break;
      }

      // iteration limit check
      if (intcount == itermax) {
        if(cxerr < error && std::abs(eerr) < error) break;
        if (std::abs(eerr) > error) {
          std::cerr<<"Error: maximum iteration step number "<<itermax<<" reached, but energy error "<<eerr<<" is larger than criterion "<<error<<std::endl;
        } else {
          std::cerr<<"Error: maximum iteration step number "<<itermax<<" reached, but phase error "<<cxerr<<" is larger than criterion "<<error<<std::endl;
        }          
        abort();
      }
      
      if (intcount>0) {
        // reset the initial data (t, B, w, X, V, Ekin)
        restore(d0);
        Ekin = Ekin0;
        // reset velocity to get correct w
        if (pars->beta>0) resolve_V();
      }

      // Dense output
      if (ip_flag) {
        // middle difference case: difference order from 1 to 2*intcount+2 (2*kappa-2; kappa=intcount+1), first one is used to storage f(x)        
        if(pars->exp_sequence==3) ndmax[intcount] = 2*intcount+3;
        // edge difference case: difference order from 1 to intcount+1, first store f(x)
        else ndmax[intcount] = intcount+2;
        
        // pd[][*]: * storage f(x) and difference
        pd[intcount] = new double*[ndmax[intcount]];
        // pd[][][*]: * storage data (t, B, w, X, V)
        for (std::size_t j=0;j<ndmax[intcount];j++) pd[intcount][j] = new double[pnn];
      }
      else pd[intcount] = NULL;

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
        // Using Romberg method
        if (method==1) EP::polynomial_extrapolation(dn,dtemp,step,dsize,intcount);
        // Using Rational interpolation method
        else EP::rational_extrapolation(dn,dtemp,step,dsize,intcount);

        // get error estimation
        double ermax=EP::extrapolation_error(dn,dsize,intcount);
        double dsfactor = EP::H_opt_factor(ermax,error,intcount);
        double werrn = (double)itercount / dsfactor;
        if (ermax>0&&werrn<werrmax) {
          werrmax = werrn;
          dsn = dsfactor;
        //dsn = std::min(dsn,0.9); // not larger than 1.0
        //dsn = std::max(dsn,pars->dtmin); // not too small
#ifdef DEBUG
          std::cerr<<"ERR factor update: iterindex="<<intcount<<"; modify factor="<<1.0/dsfactor<<"; ermax="<<ermax<<std::endl;
#endif
        }
      
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
        eerr = (Ekin-Pot+B)/B;
        std::memcpy(CX,CXN,3*sizeof(double));
      }
      
      intcount++;
#ifdef DEBUG
      std::cerr<<std::setprecision(6)<<"Iter.= "<<intcount<<" Dep.= "<<step[intcount]<<" P-err.= "<<cxerr;
      std::cerr<<" E-err.="<<(Ekin-Pot+B)/B<<" B ="<<std::setprecision(12)<<B<<std::endl;
#endif 
    }

    if (intcount+1<pars->opt_iter) dsn = std::pow(ds, (2*intcount+3)/(double)(2*pars->opt_iter+1)-1);

    // for dense output
    if (toff>0&&toff<t&&std::abs(toff-t)>pars->dterr) {

#ifdef DEBUG
      std::cerr<<"ds="<<ds<<" step[0]="<<step[0]<<" terr="<<toff-t<<" t="<<t<<" toff="<<toff<<std::endl;
#endif
      // calculate derivates from differences
      for (std::size_t i=0; i<intcount; i++) {
        // first difference pointer
        double **pdptr;
        int dpsize;
        double h;
        if(pars->exp_sequence==3) {
          // middle difference case first element is f(x)
          pdptr=&pd[i][2];
          // differece order number should be reduced by one
          dpsize = ndmax[i]-2;
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
      if(pars->exp_sequence==3) {

        npoints=3;
        
        xpoint=new double[npoints];
        xpoint[0]=0.0;
        xpoint[1]=0.5*ds;
        xpoint[2]=ds;

        nlev=new int[npoints];
        nlev[0] = 1;
        nlev[1] = ndmax[intcount-1];
        nlev[2] = 1;

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
          dfptr[i][1]=pd[intcount-1][i+1];
          dfptr[i][2]=NULL;
        }
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
      do {
        if (rf_method) {
          dsm = (dsi[0]*(tsi[1]-toff)-dsi[1]*(tsi[0]-toff))/(tsi[1]-tsi[0]);  // Use regula falsi method to find accurate ds
        }
        else {
          if(std::abs(tpre-toff)<dterr3) rf_method=true;
          dsm = (dsi[0]+dsi[1])*0.5;      // Use bisection method to get approximate region
        }
        
        EP::Hermite_interpolation_polynomial(dsm,&tpre,&pcoff,xpoint,1,npoints,nlev);
        if (tpre > toff) {
          dsi[1] = dsm;
          tsi[1] = tpre;
        }
        else {
          dsi[0] = dsm;
          tsi[0] = tpre;
        }
#ifdef DEBUG
        std::cerr<<"Find root: dsm="<<dsm<<"; t="<<tpre<<"; error="<<tpre-toff<<"; ds="<<ds<<std::endl;
#endif
      } while (std::abs(tpre-toff)>0.1*dterr);

      // Get the interpolation result
      //EP::Hermite_interpolation_polynomial(dsm,dtemp,&pcoff,xpoint,1,npoints,nlev);

      // update time factor
      dsn = -(dsm/ds);
      
      // update the results
      //restore(dtemp);

      // reset the data to original
      restore(d0);
      Ekin = Ekin0;

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
//      for (std::size_t i=0; i<=1000; i++) {
//        dsm = ds/1000*i;
//        EP::Hermite_interpolation_polynomial(dsm,&tpre,&pcoff,xpoint,1,npoints,nlev);
//        EP::Hermite_interpolation_polynomial(dsm,dtemp,pcoff,xpoint,dsize,npoints,nlev);
//        // update the results
//        restore(dtemp);
//        // resolve particle
//        resolve_XV();
//        // recalculate the energy
//        calc_Ekin();
//        // force, potential and W
//        calc_rAPW(force);
        
//        std::cerr<<"Loop: "<<dsm;
//        std::cerr<<" "<<std::setprecision(15)<<tpre;
//        std::cerr<<" "<<(Ekin-Pot+B)/B;
//        for (std::size_t i=0;i<dsize;i++) std::cerr<<std::setprecision(10)<<" "<<dtemp[i];
//        std::cerr<<std::endl;
//      }
      
#endif

      // resolve particle
      resolve_XV();
      // recalculate the energy
      // calc_Ekin();
      // force, potential and W
      calc_rAPW(force);

      // clear memory
      for (std::size_t i=0; i<intcount; i++) delete[] pn[i];
      for (std::size_t i=0; i<(std::size_t) nlev[1]-1; i++) delete dfptr[i];
      delete[] dfptr;
      delete[] pcoff;
      delete[] xpoint;
      delete[] nlev;
      delete[] fpoint;
    }

    // update chain link order
    if(num>2) update_link();

    // clear memory
    for (std::size_t i=0; i<itermax; i++) {
      delete[] dn[i];
    }
    for (std::size_t i=0; i<intcount; i++) {
      if (pd[i]!=NULL) {
        for (std::size_t j=0; j<ndmax[i]; j++) delete[] pd[i][j];
        delete[] pd[i];
      }
    }
    
#ifdef TIME_PROFILE
    profile.initstep(itermax+1);
    profile.stepcount[intcount+1]++;
    profile.t_ep += get_wtime();
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
    return -Pot;
  }

  //! Get current system energy
  /*! \return current system energy #B (negative value for bounded systems)
  */
  double getB() const {
    return B;
  }
  
  //! Get current time transformation parameter #w
  /*! #w: \f$ \frac{dw}{dt} = \sum_k \frac{\partial W}{\partial \vec{r_k}} \bullet \vec{v_k} \f$
      (see W in getW();)
       \return current time transformation parameter
  */
  double getw() const {
    return w;
  }

  //! Get current time transformation parameter #W
  /*! #W: \f$ W = \sum_{i<j} w_{ij} / |r_{ij}| (i<j) \f$ 
      - if smooth mass coefficient is used (ARC::chainpars.setM()):  \f$ w_{ij} = mm2 = \sum_{i<j} m_i  m_j /(N (N-1)/2) \f$ if \f$ m_i m_j < epi \bullet mm2 \f$ and \f$ w_{ij} = 0 \f$ otherwise.
      - else: \f$ w_{ij} = m_i m_j \f$ 

      \return current time transforamtion parameter #W
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
             <<"Time momentum B:"<<std::setw(width)<<B<<std::endl
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
template <class particle>
class chainlist{
  std::size_t num; //!< number of current particles in the list #p
  std::size_t nmax; //!< maximum number of particles allown to store
  std::size_t nchain;  //!< number of chain type members
  bool* cflag; //!< flag array to indicate whether the corresponding particle in #p is chain (true: chain; false: Particle)
  void** p; //!< particle list array (void pointer array)

public:
  //! Constructor 
  /*! Set current particle number to zero, need to use init() to allocate memory for storing particle addresses
   */
  chainlist(): num(0), nmax(0), nchain(0) {};
  
  //! Constructor with maximum number of particle \a n
  /*! Set maximum particle number to \a n and allocate memory for #p to storing the particle addresses (maximum \a n)
   */
  chainlist(const std::size_t n) { init(n); }

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

  //! Clear function
  /*! Free dynamical memory space used in particle address list #p
   */
  void clear() {
    if (nmax>0) {
      nmax = 0;
      num = 0;
      delete[] cflag;
      delete[] p;
    }
  }

  //! Destructor
  ~chainlist() {
    if (nmax>0) {
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

  //! Get number of chain members
  /*! \return Number of chain members in particle address list #p
   */
  std::size_t getNchain() const {
    return nchain;
  }

  //! Add new particle
  /*! Add new particle address at the end of the particle address list #p
    @param [in] a: new particle
   */
  void add(particle &a) {
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
    if (num+n<=nmax) {
      for (std::size_t i=0;i<n;i++) {
        cflag[num+i] = false;
        p[num+i] = &a[i];
      }
      num +=n;
    }
    else {
      std::cerr<<"Error: chainlist overflow! maximum number is "<<nmax<<", try to add "<<n<<std::endl;
      abort();
    }
  }

  //! Add a chain in particle address list #p
  /*! Add one chain's address at the end of particle address list #p
    @param [in] a: new chain particle
   */
  void add(chain<particle> &a) {
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
  particle &operator [](const std::size_t i){
    if (i>=num) {
      std::cerr<<"Error: the required index "<<i<<" exceed the current particle list boundary (total number = "<<num<<std::endl;
      abort();
    }
    if (cflag[i]) {
      return ((chain<particle>*)p[i])->cm;
    }
    else {
      return *((particle*)p[i]);
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
  chain<particle> *getSub (const std::size_t i) const {
    if (cflag[i]) return (chain<particle>*)p[i];
    else return NULL;
  }

};

}
