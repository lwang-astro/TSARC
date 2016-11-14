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

#ifdef TIME_PROFILE
#include <sys/time.h>
static double get_wtime(){
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return tv.tv_sec + 1.e-6 * tv.tv_usec;
}

struct timeprofile{
public:
  double t_apw;    // APW calculation
  double t_uplink; // update link
  double t_lf;     // leap-frog
  double t_ep;     // extrapolation
  double t_pext;   // perturber force

  int* stepcount;  // iteration step count in extrapolation

  timeprofile() {reset_tp();}

  void initstep(const std::size_t n) {
    if (!stepcount) {
      stepcount=new int[n];
      for (std::size_t i=0; i<n; i++) stepcount[i]=0;
    }
  }
  
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

/* pair AW: acceleration (potential) & dW/dr (W) function function pointer from j to i
   argument: A: acceleration vector
             P: potential
             dW: dW/dr (used for beta>0) 
             W: time transformation function (used for beta>0)
             X: relative position (1:3)
             mi: particle i mass;
             mj: particle j mass;
             mm2: smooth mass coefficient \sum m_i * m_j /(N(N-1)/2) (i<j) (Notice only calculated when m_smooth is true in chainpars)
             epi: adjustable parameter (set in chainpars)
 */
typedef void (*pair_AW) (double *, double &, double *, double &, const double *, const double &, const double &, const double &, const double &);

/* pair Ap: acceleration calculation function pointer from j to i
   argument: A: acceleration vector
             P: potential
             pi: position vector i
             pj: position vector j
             mi: particle mass i;
             mj: particle mass j;
 */
typedef void (*pair_Ap) (double *, double &, const double*, const double*, const double&, const double&);

/* Newtonian AW from k to j
   mm2 = sum (i<j) m_i*m_j/(N(N-1)/2)
*/
void Newtonian_AW (double A[3], double &P, double dW[3], double &W, const double xjk[3], const double &mj, const double &mk, const double &mm2, const double &epi) {

  //distance
  double rjk = std::sqrt(xjk[0]*xjk[0]+xjk[1]*xjk[1]+xjk[2]*xjk[2]);

  // mass parameters
  double mmjk = mj*mk;
  double wjk;
  if (mm2>0 && epi>0) {
    // Wjk = mm2 if m_i*m_j < epi*m'^2; 0 otherwise;
    if (mmjk<epi*mm2) wjk = mm2;
    else wjk = 0;
  }
  else {
    // Wjk = m_i*m_j
    wjk = mmjk;
  }
  
  //Potential energy==================================//
  P = mmjk / rjk;
  //Transformation coefficient========================//
  W = wjk / rjk;
        
  //Acceleration======================================//
  double rjk3 = rjk*rjk*rjk;
  double mor3 = mk / rjk3;
  A[0] = mor3 * xjk[0];
  A[1] = mor3 * xjk[1];
  A[2] = mor3 * xjk[2];

  //d W / d r=========================================//
  mor3 = wjk / rjk3;
  dW[0] = mor3 * xjk[0];
  dW[1] = mor3 * xjk[1];
  dW[2] = mor3 * xjk[2];
  
}

/* Newtonian force from p to i*/
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
 
// chain parameter controller

class chainpars{
template <class T> friend class chain;
private:
  // pair force
  pair_AW pp_AW;
  pair_Ap pp_Ap;
  
  //time step integration parameter
  double alpha; // logH cofficient
  double beta;  // TTL cofficient
  double gamma; // constant

  //  mass coefficients parameter
  double m_epi; // smooth parameter
  bool m_smooth; // whether to use smooth mass coefficients

  // time step
  double dtmin; // minimum physical time step
  double dterr; // physical time error criterion

  // extrapolation control parameter
  double exp_error; // relative error requirement for extrapolation
  std::size_t exp_itermax; // maximum times for iteration.
  int exp_methods;  // 1: Romberg method; others: Rational interpolation method
  int exp_sequences; // 1: even sequence {h, h/2, h/4, h/8 ...}; 2: Bulirsch & Stoer sequence {h, h/2, h/3, h/4, h/6, h/8 ...}; other. 4k sequence {h/2, h/6, h/10, h/14 ...}

  int* step; //substep sequence
  std::size_t  opt_iter; // optimaized iteration index

  int** bin_index; // binomial coefficients

public:

  chainpars(): alpha(1.0), beta(0.0), gamma(0.0), m_epi(0.001), m_smooth(true) {
    step = NULL;
    bin_index = NULL;
    setEXP(1E-10, 5.4E-20, 1E-6, 20, 2, 2, 5);
    pp_AW = &Newtonian_AW;
    pp_Ap = &Newtonian_Ap;
  }

  chainpars(pair_AW aw, pair_Ap ap, const double a, const double b, const double g, const double e=0.001, const bool mm=true, const double error=1E-10, const double dtm=5.4e-20, const double dte=1e-6, const std::size_t itermax=20, const int methods=2, const int sequences=2, const std::size_t optiter=5) {
    step = NULL;
    bin_index = NULL;
    setabg(a,b,g);
    setM(e,mm);
    setEXP(error,dtm,dte,itermax,methods,sequences,optiter);
    setA(aw,ap);
  }

  ~chainpars() {
    if (step!=NULL) delete[] step;
    if (bin_index!=NULL) {
      for (std::size_t i=0;i<exp_itermax;i++) 
        if (bin_index[i]!=NULL)
          delete[] bin_index[i];
      delete[] bin_index;
    }
  }

  /* setA
     function: set pair acceleration function
     argument: aw: pair_AW
               ap: pair_Ap
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
  
     
  /* setabg
     function: set time step integration parameter alpha, beta, gamma
     argument: a: alpha (logH)
               b: beta  (TTL)
               g: gamma (constant)
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

  /* setM
     function: set mass coefficient parameters epi and smooth_flag
     argument: e: epi   (smooth parameter)
               option: m_smooth (whether use smooth mass coefficients)
  */
  void setM(const double e=0.001, const bool mm=true) {
    m_epi = e;
    m_smooth = mm;
    if (m_epi==0&&m_smooth) {
      std::cerr<<"Error: smooth mass coefficients are used, but smooth coefficient epi is zero!\n";
      abort();
    }
    if (alpha==0&&gamma==0) {
      std::cerr<<"Warning: alpha=0 and gamma=0, the smooth mass coefficients may cause initial w = 0 and result in zero time step!";
    }
  }

  /* setEXP
     function: set extrapolation parameters
     argument: error: relative error requirement for extrapolation
               dtmin: minimum physical time step
               itermax: maximum times for iteration.
               methods: 1: Romberg method; others: Rational interpolation method
               sequences: 1: even sequence {h, h/2, h/4, h/8 ...}; 2: Bulirsch & Stoer sequence {h, h/2, h/3, h/4, h/6, h/8 ...}; 3: 4k sequence {h, h/2, h/6, h/10, h/14 ...}; others: Harmonic sequence {h, h/2, h/3, h/4 ...}
  */
  void setEXP(const double error=1E-10, const double dtm=5.4e-20, const double dte=1e-6, const std::size_t itermax=20, const int methods=2, const int sequences=2, const std::size_t optiter=5) {
    exp_error = error;
    exp_methods = methods;
    exp_sequences = sequences;
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
  
  /*  const double &alpha()  const { return S_alpha;  }
  const double &beta()   const { return S_beta;   }
  const double &gamma()  const { return S_gamma;  }
  const double &m_epi()    const { return M_epi;    } 
  const double &m_smooth() const { return M_smooth; }
  const double &exp_error()  const { return EXP_error;}
  const std::size_t &exp_itermax() const { return EXP_itermax; }
  const int &exp_methods()   const { return EXP_methods; }
  const int &sequences()     const { return EXP_sequences; }
  const double &dtmin()      const { return INT_dtmin; }*/
  
  
};

// Chain class
template <class particle>
class chain{
  typedef double double3[3];
  double3 *X;  // relative position
  double3 *V;  // relative velocity
  std::size_t *list;   // chain index list
  double3 *acc; // acceleration
  double3 *pf;  // perturber force
  double3 *dWdr; // \partial Omega/ \partial rk

  //integration parameters=======================================//
  double t;     //time
  double w;    //time transformation parameter
  double B;    //Binding energy (time momentum)

  //template parameters==========================================//
  double W;     //time transformation function
  double Ekin;  //kinetic energy
  double Pot;   //potential
  double mm2;  // mean mass production \sum m_i*m_j/(N(N-1)/2) (i<j)

  //number =======================================================//
  std::size_t num;      //total number of chain particles
  std::size_t nmax;     //maximum number 

  //monitor flags
  bool F_Pmod;     // indicate whether particle list is modified (true: modified)
  int  F_Porigin; // indicate whether particle is shifted back to original frame (1: original frame: 0: center-of-mass frame; 2: only position is original frame)

  // chain parameter controller
  const chainpars *pars;

public:

  //center mass=======================================//
  particle cm;

  //  double cm.getMass();  //total mass
  //  double cm.getPos();   // position center
  //  double cm.getVel();   // velocity center

  //particle list=====================================//
  chainlist<particle> p;

  //perturber list====================================//
  chainlist<particle> pext;

#ifdef TIME_PROFILE
  timeprofile profile;
#endif
  
  //initialization
  //  template <class particle>
  chain(std::size_t n, const chainpars &par):  pars(&par) {
    nmax=0;
    mm2=-std::numeric_limits<double>::infinity();
    allocate(n);
  }

  chain(const chainpars &par): pars(&par), F_Pmod(false), F_Porigin(1), num(0), nmax(0) {
    mm2=-std::numeric_limits<double>::infinity();
  }

  // re-allocate function
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

  // clear function
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

  //destruction
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
  /* update_num
     function: to update the num due to current particle list number
     argument: n: particle number
  */
  void update_num(const std::size_t n) {
    if (n>nmax) {
      std::cerr<<"Error: particle number "<<n<<" is larger than Chain number limit "<<num<<std::endl;
      abort();
    }
    else{
      // Update Chain number
      num = n;
    }
  }

  /* generate_list  ============================
     function: generate chain list by N^2 searching particle lists
  */
  //  template <class particle>
  void generate_list() {
    bool *is_checked=new bool[num];
    for (std::size_t i=0; i<num; i++) is_checked[i] = false;
    std::size_t inext=0;
    for (std::size_t i=0; i<num; i++) {
      //initial rjk=======================================//
      //mark checked particle=============================//
      is_checked[inext] = true;

      //initial chain_mem=================================//
      list[i]=inext;
      std::size_t inow=inext;
    
      //make chain========================================//
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

  /* calc_XV
     function: get chain member relative position X and velocity V based on current list
  */
  //  template <class particle>
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


  /* calc_mm2
     function: Get averaged mass coefficients for calculating Wjk
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

  
  /* calc_rAPW
     function: Get distance Matrix, acceleration, dm/dr, potential and transformation parameter
               based on particle mass (obtain from p), current X & V and Wjk
               (notice the acceleration array index k and the distance matrix index j,k order follow particle p to avoid additional shift when chain list change)
     argument: force: external force (acceleration) for each particle, (not perturber forces)
               resolve_flag: flag to determine whether to resolve sub-chain particles for force calculations. 
  */
  void calc_rAPW (const double3 *force=NULL, const bool resolve_flag=false) {
#ifdef TIME_PROFILE
    profile.t_apw -= get_wtime();
#endif
    // template xjk vector
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
        //Acceleration==========================================//
        acc [lj][k]=0.0;
        //dw/dr ================================//
        dWdr[lj][k]=0.0;
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

            //Acceleration======================================//
            At[0] += Atemp[0];
            At[1] += Atemp[1];
            At[2] += Atemp[2];
            //Potential=========================================//
            if (k>j) Pt += Ptemp;
          }
          // center shift back
          ck->center_shift_X();
        }

        //Acceleration======================================//
        acc[lj][0] += At[0];
        acc[lj][1] += At[1];
        acc[lj][2] += At[2];

        //d W / d r=========================================//
        dWdr[lj][0] += dWt[0];
        dWdr[lj][1] += dWt[1];
        dWdr[lj][2] += dWt[2];

        if (k>j) {
          //Potential energy==================================//
          Pot_c += Pt;
          //Transformation coefficient========================//
          W_c += Wt;
        }
      }
      //add external acceleration======================================//
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

  /* calc_Ekin
     function: calculate kinetic energy
  */
  //  template <class particle>
  void calc_Ekin(){
    Ekin = 0.0;
    for (std::size_t i=0; i<num; i++) {
      const double *vi=p[i].getVel();
      Ekin += 0.5 * p[i].getMass() * (vi[0]*vi[0]+vi[1]*vi[1]+vi[2]*vi[2]);
    }
  }

  /* calc_dt_X
     function: calculate physical time step dt for X based on ds
     argument:  ds:  step size s (not physical time step) 
     return:    dt:  physical integration time step for X
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
  
  /* calc_dt_V
     function: calculate physical time step dt for V based on ds
     argument:  ds:  step size s (not physical time step) 
     return:    dt:  physical integration time step for V
  */
  double calc_dt_V(const double ds) {
    // determine velocity integration time step
    return ds / (pars->alpha * Pot + pars->beta * W + pars->gamma);
  }
  
  /* step_forward_X
     function: one step integration of X 
     argument: dt: physical time step dt for X
   */
  void step_forward_X(const double dt) {
    // step forward relative X
    for (std::size_t i=0;i<num-1;i++) {
      X[i][0] += dt * V[i][0];
      X[i][1] += dt * V[i][1];
      X[i][2] += dt * V[i][2];
    }
  }

  /* step_forward_V
     function: one step integration of V
     argument: dt: physical time step dt for V
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

  /* step_forward_Bw
     functrion: one step integration of B and w.
                B += dt * \sum ( - m_k * <v_k> \dot f_k);
                w += dt * \sum ( dm/dr_k \dot v_k);
     argument: dt: time step for V
               ave_v: averaged velocity
               force: external force
               p: particle list (only use mass)
               fpf: if perturber force is not zero, it is true
   */
  //  template <class particle>
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

  /* resolve_XV
     function: resolve relative X, V to physical x, v and calculated the averaged velocity of old and new values.
               Notice the center-of-mass particle mass in Chain.cm is used.
               The total mass of particles should be consistent with cm.getMass(). Otherwise update Chain.cm first.
     argument: ave_v: averaged velocity array (return values)
   */
  //  template <class particle>
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

      //center-of-mass position and velocity
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

  /* resolve_X
     function: resolve relative X to physical x (center-of-mass frame)
               Notice the center-of-mass particle mass in Chain.cm is used.
               The total mass of particles should be consistent with cm.getMass(). Otherwise update Chain.cm first.
   */
  //  template <class particle>
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

      //center-of-mass position and velocity
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

  /* resolve_V
     function: resolve relative V to physical  v (center-of-mass frame)
               Notice the center-of-mass particle mass in Chain.cm is used.
               The total mass of particles should be consistent with cm.getMass(). Otherwise update Chain.cm first.
   */
  //  template <class particle>
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

      //center-of-mass position and velocity
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

  /* pert_force
    function: get perturber force
    argument: resolve_flag: whether resolve perturber chain
    return:   true: pertubers exist. false: no perturbers
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

              //Acceleration======================================//
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
     
  /* update_link
     function: update chain link based on current relative distance matrix and V
     return:   if link is modified, return true
   */
  bool update_link(){
#ifdef TIME_PROFILE
    profile.t_uplink -= get_wtime();
#endif
    // indicator
    bool modified=false;
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

  /* center_shift ==================================
     function: shift positions and velocities of N (num) particles (p) based on their center-of-mass, write center-of-mass particle to chain
     p: particle list (read data and return shifted list);
  */
  void center_shift_init() {
    //center mass=======================================//
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
    cm.setPos(cmr);
    cm.setVel(cmv);

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

  /* center_shift_inverse_X==================================
     function: shift back positions of N (num) particles (p) based on chain center-of-mass
     p: particle list (read data and return shifted list);
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

  /* center_shift_X==================================
     function: shift positions and velocities of N (num) particles (p) based on chain center-of-mass
     p: particle list (read data and return shifted list);
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


  /* mid_diff_calc
    function: central difference cumulative calculator for t only
              If i is 0, the value will be reset to zero, and for the same n, this function need to be called n+1 times to complete the n'th order difference
    argument: tdiff: one dimensional array to storage the difference of t for different order. Array size is [n], If tdiff is NULL, no calculation will be done.
              n: maximum difference level (count from 1)
              i: current position of t_i corresponding to the binomial coefficients (count from 0 as starting point)
              ndiv: current substep number can be used
              tflag: calculate time only if true
  */
  void mid_diff_calc(double **dpoly, const std::size_t nmax, const std::size_t i, const int ndiv, bool tflag=true) {
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
        if (!tflag) {
          const std::size_t dsize=6*num-3;
          //initial dpoly to zero
          for (std::size_t j=0; j<nmax; j++) if(i==0) std::memset(dpoly[j],0,dsize*sizeof(double));
        }
        else if(i==0) for (std::size_t j=0; j<nmax; j++) dpoly[j][0] = 0.0;
        
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
              dpoly[j][0] += coff * t;
              if (!tflag) {
                dpoly[j][1] += coff * B;
                dpoly[j][2] += coff * w;
                for (std::size_t k=0; k<num-1; k++) {
                  for (std::size_t kk=0; kk<3; kk++) {
                    dpoly[j][3*(1+k)+kk] += coff * X[k][kk];
                    dpoly[j][3*(num+k)+kk] += coff * V[k][kk];
                  }
                }
              }
            }
          }
        }
      }
    }
  }


  /* edge_diff_calc
     function: left forward and right backward difference cumulative calculator for t, B, w, X, V
              If i is 0, the array will be reset to zero, and for the same n, this function need to be called n+1 times (i=0...n+1) to complete the n'th order difference.
    argument: dpoly: two dimensional array to storage the results. Array size is [n][12*num-6+1] (last is for safety checking). first [] indicate different level of differences, second [] indicate the left and right difference of data (t, B, w, X, V), similar as backup/restore data. If dpoly is NULL, no calculation will be done.    
              nmax: maximum difference level (count from 1)
              i: current position of data (t, B, w, X, V) corresponding to the binomial coefficients (count from 0 as starting poin              ndiv: total step number
              tflag: calculate time only if true
t)
   */
  void edge_diff_calc(double **dpoly, const int nmax, const std::size_t i, const int ndiv, bool tflag=true) {
    // safety check
    if (dpoly!=NULL) {
      // safety check
      if (!tflag) {
        const std::size_t dsize=12*num-6;
        // i indicate the position, i=0 means x0
        if(i==0) for (int j=0; j<nmax; j++) std::memset(dpoly[j],0,dsize*sizeof(double)); // reset data array at i=0
      }
      else if(i==0) for (int j=0; j<nmax; j++) {
          dpoly[j][0]=0.0;
          dpoly[j][1]=0.0;
        }
      
      if (nmax>ndiv) {
        std::cerr<<"Error: maximum difference order "<<nmax<<" > total step number "<<ndiv<<"!\n";
        abort();
      }

      if (i<=nmax||i>=ndiv-nmax) {
        int** binI = pars->bin_index;
        for (int j=0; j<nmax; j++) {
          // j+1 indicate the difference degree, count from 1 (first difference)
          const int n = j+1;
              
          // left edge, forward difference: (-1)^(n-i) (n n-i) f(x0+i*h) (i from 0 to n)
          int ileft = (int)(n-i);   // n-i
          if (ileft>=0) {
            double coff = ((ileft%2)?-1:1)*binI[n][ileft];
            dpoly[j][0] += coff * t;
            if (!tflag) {
              dpoly[j][1] += coff * B;
              dpoly[j][2] += coff * w;
              for (std::size_t k=0; k<num-1; k++) {
                for (std::size_t kk=0; kk<3; kk++) {
                  dpoly[j][3*(1+k)+kk] += coff * X[k][kk];
                  dpoly[j][3*(num+k)+kk] += coff * V[k][kk];
                }
              }
            }
#ifdef DEBUG
            std::cerr<<"Poly left: n="<<n<<" ik="<<ileft<<" i="<<i<<" coff="<<coff<<" t="<<t<<std::endl;
#endif
          }

          // right edge, backward difference: (-1)^(n-i) (n n-i) f(xn-(n-i)*h) (i from 0 to n)
          int ishift = ndiv-n;
          int iright= ileft+ishift;     // n-i
          if (i>=ishift) {
            double coff = ((iright%2)?-1:1)*binI[n][iright];
            if (!tflag) {
              const std::size_t ir = 6*num-3;
              dpoly[j][0+ir] += coff * t;
              dpoly[j][1+ir] += coff * B;
              dpoly[j][2+ir] += coff * w;
              for (std::size_t k=0; k<num-1; k++) {
                for (std::size_t kk=0; kk<3; kk++) {
                  dpoly[j][3*(1+k)+kk+ir] += coff * X[k][kk];
                  dpoly[j][3*(num+k)+kk+ir] += coff * V[k][kk];
                }
              }
            } else dpoly[j][1] += coff * t;
#ifdef DEBUG
            std::cerr<<"                 Poly right: n="<<n<<" ik="<<iright<<" i="<<i<<" coff="<<coff<<std::endl;
#endif
          }
        }
      }
    }
  }

  /* edge_dev_calc
     function: calcualte derivates from differences: \d^(n)f/ h^n
     argument: dpoly: two dimensional storing differences, will be updated to derivates. Array size is [n][12*num-6+1] (last is for safety checking). first [] indicate different level of differences, second [] indicate the left and right difference of data (t, B, w, X, V), similar as backup/restore data. If dpoly is NULL, no calculation will be done.
               h: step size
               nmax: maximum difference order
               dsize: data number
   */
  void edge_dev_calc(double **dpoly, const double h, const int nmax, const int dsize) {
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
  /* add particle
     function: add particle into chainlist p
     argument: particle a
   */
  void addP(particle &a) {
    if (F_Porigin!=1) std::cerr<<"Warning!: particle list are (partically) in the center-of-mass frame, dangerous to add new particles!\n";
    p.add(a);
    F_Pmod=true;
  }

  void addP(chain<particle> &a) {
    if (F_Porigin!=1) std::cerr<<"Warning!: particle list are (partically) in the center-of-mass frame, dangerous to add new particles!\n";
    p.add(a);
    F_Pmod=true;
  }
  
  void addP(const std::size_t n, particle a[]) {
    if (F_Porigin!=1) std::cerr<<"Warning!: particle list are (partically) in the center-of-mass frame, dangerous to add new particles!\n";
    p.add(n,a);
    F_Pmod=true;
  }

  /* remove particle
     function: remove particle from p
     argument: i: particle index
               option: true: shift last particle to current position (defaulted);
                       false: shift all right particle to left by one
  */
  void removeP(const std::size_t i, bool option=true) { p.remove(i,option); F_Pmod=true; }


  /* initPext
     function: allocate array for Pext list
     argument: n: number of perturbers
   */
  void initPext(const std::size_t n) {
    if (pext.getN()) {
      std::cerr<<"Error: Perturber list is already initialized!\n";
      abort();
    }
    pext.init(n);
  }

  
  /* add perturber particle
     function: add pertuber particle into chainlist p
     argument: particle a
   */
  void addPext(particle &a) { pext.add(a);}
  void addPext(chain<particle> &a) { pext.add(a);}
  void addPext(const std::size_t n, particle a[]) { pext.add(n,a); }

  /* remove particle
     function: remove particle from p
     argument: i: particle index
               option: true: shift last particle to current position (defaulted);
                       false: shift all right particle to left by one
  */
  void removePext(const std::size_t i, bool option=true) { pext.remove(i,option); }
  
  /* isPmod
     function: reture true if particle list is modifed, in this case, the chain may need to be initialized again
     return: true/false
  */
  bool isPmod() const { return F_Pmod; }

  /* isPorigin
     function: return true if particle list are in original frame
     return: 1 (original frame),2 (position in original frame only) ,0 (center-of-mass frame)
  */
  int isPorigin() const { return F_Porigin; }
  
  /* init
     function: initialization chain from particle p
     argument: time: current time
     force: external force
  */
  void init(const double time, const double3* force=NULL) {
    //update number indicator
    update_num(p.getN());
    
    //Generate chain link list
    generate_list();

    // calculate perturber forces
    //    pert_force();
    
    //set center-of-mass
    if (F_Porigin==1) {
      center_shift_init();
      F_Porigin=0;
    }
    else {
      std::cerr<<"Error: particles are not in original frame!\n";
      abort();
    }

    //set member relative position and velocity
    calc_XV();
    
    // if smooth mass coefficients are used, calculate mm2;
    if (pars->m_smooth) calc_mm2();

    //set relative distance matrix, acceleration, potential and transformation parameter
    calc_rAPW(force);

    //Initial intgrt value t
    t = time;

    //kinetic energy
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


  /* center_shift_inverse==================================
     function: shift back to original positions and velocities of N (num) particles (p) based on their center-of-mass
     p: particle list (read data and return shifted list);
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

  /* center_shift==================================
     function: shift positions and velocities of N (num) particles (p) based on their center-of-mass
     p: particle list (read data and return shifted list);
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

  /* backup
     function: backup integration data to one dimensional array
     argument: db: backup array (size should be 6*num-3)
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
     
  /* restore
     function: restore integration data from one dimensional array
     argument: db: data array (size should be 6*num-3)
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

  /* Leapfrog_step_forward
     function: integration n steps, particles' position and velocity will be updated in the center-of-mass frame
     argument: s: whole step size
               n: number of step division, real step size = (s/n), for Leapfrog integration, it is X(s/2n)V(s/n)X(s/n)V(s/n)..X(s/2n)
               force: external force (not perturber forces which are calculated in pert_force)
               check_flag: 2: check link every step; 1: check link at then end; 0 :no check
               dpoly: interpolation polynomial coefficients array (size shold be [n][(2*nrel+1)*3])
               ndmax: maximum difference limit
///               recur_flag: flag to determine whether to resolve sub-chain particles for force calculations. notice this require the sub-chain to be integrated to current physical time. Thus is this a recursion call (tree-recusion-integration)
///               upforce: void (const particle * p, const particle *pext, double3* force). function to calculate force based on p and pext, return to force
  */             
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
    bool fpf = false; // perturber force indicator
    const int np = pext.getN();
    if (np>0) fpf = true;

    //for polynomial coefficient calculation first point
    // middle difference (first array is used to store the f_1/2)
    if(pars->exp_sequences==3) mid_diff_calc(&dpoly[1],ndmax-1,0,n);
    // edge difference
    else edge_diff_calc(dpoly,ndmax,0,n); 
                                               
    
    //integration loop
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
      if(pars->exp_sequences==3) {
        if (i==n/2) dpoly[0][0]=t; // f(x)
        mid_diff_calc(&dpoly[1],ndmax-1,i+1,n);
      }
      // edge difference
      else edge_diff_calc(dpoly,ndmax,i+1,n);
    }

    // resolve X at last, update p.x (dependence: X)
    resolve_X();

    // Update rjk, A, Pot, dWdr, W (notice A will be incorrect since pf is not updated)
    calc_rAPW(force);

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

  /* extrapolation_integration
     function: extrapolation method to get accurate integration
     argument: ds: whole step size
               toff: ending physical time (negative means no ending time and integration finishing after n steps)
               force: external force (not perturber forces which are calculated in pert_force)
     return:   new step size modification factor
               if it is negative and toff>0; the abs of value is factor to be used for redoing the extrapolation_integrration of this step to reach the toff. And also in this case the data will not be integrated.
   */
  double extrapolation_integration(const double ds, const double toff=-1.0, const double3* force=NULL) {
    // get parameters
    const double error = pars->exp_error;
    const std::size_t itermax = pars->exp_itermax;
    const int methods = pars->exp_methods;
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
    bool ip_flag = true; // interpolation coefficient calculation flag
    double** pd[itermax]; // central difference, [*] indicate different accuracy level 
    int ndmax[itermax]; //maximum difference order
    std::size_t pnn;  // data size
    if(pars->exp_sequences==3) pnn = 1;     // middle difference case
    else pnn = 2; // edge two points case

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
    
    std::size_t intcount = 0; // iteration counter
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
        // middle difference case: difference order from 1 to 2*intcount+1 (2*kappa-1; kappa=intcount+1), first one is used to storage f(x)        
        if(pars->exp_sequences==3) ndmax[intcount] = 2*intcount+2;
        // edge difference case: difference order from 1 to intcount+1
        else ndmax[intcount] = intcount+1;
        
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
        if (methods==1) EP::polynomial_extrapolation(dn,dtemp,step,dsize,intcount);
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
    if (toff>0&&toff<t) {

#ifdef DEBUG
      std::cerr<<"ds="<<ds<<" step[0]="<<step[0]<<std::endl;
#endif
      // calculate derivates from differences
      for (std::size_t i=0; i<intcount; i++) {
        // first difference pointer
        double **pdptr;
        int dpsize;
        double h;
        if(pars->exp_sequences==3) {
          // middle difference case first element is f(x)
          pdptr=&pd[i][1];
          // differece order number should be reduced by one
          dpsize = ndmax[i]-1;
          // step size
          h = 2*ds/(double)step[i];
        }
        else {
          pdptr=pd[i];
          dpsize = ndmax[i];
          h = ds/(double)step[i];
        }

        edge_dev_calc(pdptr,h,dpsize,pnn); 
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
          if (methods==1) EP::polynomial_extrapolation(pn,pd[j][i],&step[istart],pnn,j-istart);
          if (methods==2) EP::rational_extrapolation(pn,pd[j][i],&step[istart],pnn,j-istart);
#ifdef DEBUG
          std::cerr<<"Poly extra order="<<j-istart<<" step("<<istart<<")="<<step[istart]<<" t X11^("<<i+1<<")_"<<j<<"="<<pd[j][i][0]<<"\t"<<pd[j][i][3]<<std::endl;
#endif
        }
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
      if(pars->exp_sequences==3) {

        npoints=3;
        
        xpoint=new double[npoints];
        xpoint[0]=0.0;
        xpoint[1]=0.5*ds;
        xpoint[2]=ds;

        nlev=new int[npoints];
        nlev[0] = 1;
        nlev[1] = 2*(int)intcount;
        nlev[2] = 1;

        fpoint=new double*[npoints];
        fpoint[0] = d0;
        fpoint[1] = pd[intcount-1][0];
        fpoint[2] = dn[intcount-1];

        // \sum nlev = 2*intcount+2;
        pcoff = new double[2*intcount+2];

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
        nlev[0] = (int)intcount+1;
        nlev[1] = nlev[0];
        
        // store f(x)
        fpoint=new double*[npoints];
        fpoint[0] = d0;
        fpoint[1] = dn[intcount-1];

        // \sum nlev = 2*intcount+2;
        pcoff= new double[2*intcount+2];

        dfptr=new double**[intcount];
        for (int i=0;i<intcount;i++) {
          dfptr[i] = new double*[2];
          dfptr[i][0] = pd[intcount-1][i];
          dfptr[i][1] = &pd[intcount-1][i][1];
        }
      }

      
      // Hermite interpolation
      EP::Hermite_interpolation_coefficients(&pcoff,xpoint,fpoint,dfptr,1,npoints,nlev);

      // Iteration to get correct physical time position
      double dsi[2]   = {0,ds};    // edges for iteration
      double tsi[2]   = {d0[0],t}; // edges values
      double dsm,tpre;    // expected ds and t;
      const double dterr = pars->dterr;
      const double dterr3 = 1000*dterr;  //1000 * time error criterion

      bool rf_method=false;
      do {
        if (rf_method) dsm = (dsi[0]*(tsi[1]-toff)-dsi[1]*(tsi[0]-toff))/(tsi[1]-tsi[0]);  // Use regula falsi method to find accurate ds
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
      } while (std::abs(tpre-toff)>dterr);

      // Get the interpolation result
      //EP::Hermite_interpolation_polynomial(dsm,dtemp,&pcoff,xpoint,1,npoints,nlev);

      // update time factor
      dsn = -(dsm/ds);
      
      // update the results
      //restore(dtemp);
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
      /*      for (std::size_t i=0; i<=1000; i++) {
        dsm = ds/1000*i;
        EP::Hermite_interpolation_polynomial(dsm,dtemp,pcoff,xpoint,dsize,npoints,nlev);
        // update the results
        restore(dtemp);
        // resolve particle
        resolve_XV();
        // recalculate the energy
        calc_Ekin();
        // force, potential and W
        calc_rAPW(force);
        
        std::cerr<<"Loop: "<<dsm;
        std::cerr<<" "<<(Ekin-Pot+B)/B;
        for (std::size_t i=0;i<dsize;i++) std::cerr<<std::setprecision(10)<<" "<<dtemp[i];
        std::cerr<<std::endl;
        }*/
      
#endif

      // resolve particle
      resolve_XV();
      // recalculate the energy
      // calc_Ekin();
      // force, potential and W
      calc_rAPW(force);

      for (std::size_t i=0; i<intcount; i++) delete[] pn[i];
      //for (std::size_t i=0; i<dsize; i++) delete[] pcoff[i];
      delete[] pcoff;
    }

    // update chain link order
    if(num>2) update_link();

    for (std::size_t i=0; i<itermax; i++) {
      delete[] dn[i];
    }
    for (std::size_t i=0; i<intcount; i++) {
      if (pd[i]!=NULL) {
        for (std::size_t j=0; j<=i; j++) delete[] pd[i][j];
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

  /* gettime
     function: return physical time
     return: t
  */
  double getTime() const {
    return t;
  }
  /* get_Ekin
     function: get kinetic energy
     return: Ekin
  */
  double getEkin() const {
    return Ekin;
  }

  /* get_Pot
     function: get potetnial energy (negative value)
     return: Pot
  */
  double getPot() const {
    return -Pot;
  }

  /* get_B
     function: get potetnial energy (negative value)
     return: Pot
  */
  double getB() const {
    return B;
  }
  
  /* get_w
     function: get potetnial energy (negative value)
     return: Pot
  */
  double getw() const {
    return w;
  }

  /* get_W
     function: get potetnial energy (negative value)
     return: Pot
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
   */
  void set_X(const std::size_t i, const std::size_t k, const double value) {
    X[i][k] = value;
  }


  //  template <class particle>
  void update_rAPW(const double3* force) {
    calc_rAPW(force);
  }
  
  /* print
     function: print all data
     argument: width of each output value
  */
  //  template <class particle>
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

// Generalized list to store chain and Particle members
template <class particle>
class chainlist{
  std::size_t num;
  std::size_t nmax;
  std::size_t nchain;  // number of sub-chain
  bool* cflag; // flag to indicater whether it is chain (true: chain; false: Particle)
  void** p;

public:
  // initialization
  chainlist(): num(0), nmax(0), nchain(0) {};
  
  chainlist(const std::size_t n) { init(n); }

  void init(const std::size_t n) {
    num = 0;
    nmax = n;
    nchain = 0;
    cflag=new bool[n];
    p=new void*[n];
  }

  void clear() {
    if (nmax>0) {
      nmax = 0;
      num = 0;
      delete[] cflag;
      delete[] p;
    }
  }

  ~chainlist() {
    if (nmax>0) {
      delete[] cflag;
      delete[] p;
    }
  }

  // get particle number
  std::size_t getN() const {
    return num;
  }

  // get chain number
  std::size_t getNchain() const {
    return nchain;
  }

  // add new particle
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

  // add particle list
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

  // add new chain particle
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

  /* remove particle
     argument: option: true: shift last particle to current position (defaulted);
                       false: shift all right particle to left by one
               i: particle index
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

  // return the particle reference or center-of-mass particle reference if it is a chain
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

  /* ischain
     function: check whether a member is chain
     argument: i: index
     return:  true: chain; false: particle
   */
  bool isChain (const std::size_t i) const {
    return cflag[i];
  }

  chain<particle> *getSub (const std::size_t i) const {
    if (cflag[i]) return (chain<particle>*)p[i];
    else return NULL;
  }

};

