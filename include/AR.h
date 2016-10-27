#pragma once

#include <iostream>
#include <string.h>
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

  void initstep(const std::size_t n) { if (!stepcount) stepcount=new int[n];}
  
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

  std::size_t* step; //substep sequence
  std::size_t  opt_iter; // optimaized iteration index

public:

  chainpars(): alpha(1.0), beta(0.0), gamma(0.0), m_epi(0.001), m_smooth(true) {
    step = NULL;
    setEXP(1E-10, 5.4E-20, 1E-6, 20, 2, 2, 5);
    pp_AW = &Newtonian_AW;
    pp_Ap = &Newtonian_Ap;
  }

  chainpars(pair_AW aw, pair_Ap ap, const double a, const double b, const double g, const double e=0.001, const bool mm=true, const double error=1E-10, const double dtm=5.4e-20, const double dte=1e-6, const std::size_t itermax=20, const int methods=2, const int sequences=2, const std::size_t optiter=5) {
    step = NULL;
    setabg(a,b,g);
    setM(e,mm);
    setEXP(error,dtm,dte,itermax,methods,sequences,optiter);
    setA(aw,ap);
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
               sequences: 1: even sequence {h, h/2, h/4, h/8 ...}; 2: Bulirsch & Stoer sequence {h, h/2, h/3, h/4, h/6, h/8 ...}; other: 4k sequence {h, h/2, h/6, h/10, h/14 ...}
  */
  void setEXP(const double error=1E-10, const double dtm=5.4e-20, const double dte=1e-6, const std::size_t itermax=20, const int methods=2, const int sequences=2, const std::size_t optiter=5) {
    exp_error = error;
    exp_methods = methods;
    exp_sequences = sequences;
    dterr = dte;
    dtmin = dtm;
    exp_itermax = itermax;
    opt_iter=optiter;
    // reset step array
    if (step!=NULL) delete[] step;
    step = new std::size_t[itermax+1];
    // calculate sequences of steps
    // Romberg (even) sequence {h, h/2, h/4, h/8 ...}
    if (sequences==1) {
      step[0] = 1;
      for (std::size_t i=1;i<=itermax;i++) {
        step[i] = 2*step[i-1];
      }
    }
    // Bulirsch & Stoer sequence {h, h/2, h/3, h/4, h/6, h/8 ...}
    else if (sequences==2) {
      step[0] = 1;
      std::size_t stepeven = 2; 
      std::size_t stepodd = 3;
      for (std::size_t i=1;i<=itermax;i++) {
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
    // E. Hairer (4k) sequences {h, h/2, h/6, h/10, h/14 ...}
    else {
      step[0] = 1;
      step[1] = 2;
      for (std::size_t i=2;i<=itermax;i++) {
        step[i] = step[i-1] + 4;
      }
    }
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

  //integration paramenters=======================================//
  double Ekin;  //kinetic energy
  double Pot;   //potential
  double W;     //time transformation function
  double t;     //time
  double w;    //time transformation parameter
  double B;    //Binding energy (time momentum)
  double mm2;  // mean mass production \sum m_i*m_j/(N(N-1)/2) (i<j)

  //number =======================================================//
  std::size_t num;      //total number of chain particles
  std::size_t nmax;     //maximum number 

  //monitor flags
  bool F_Pmod;     // indicate whether particle list is modified (true: modified)
  int  F_Porigin; // indicate whether particle is shifted back to original frame (1: original frame: 0: center-of-mass frame; 2: only position is original frame)
  bool F_mm2s;     // indicate whether mm2 is calculated.

  //error control
  std::size_t INT_q; // last step corresponding index
  double INT_s;      // last step size
  double INT_errmax; // last step errmax

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
  chain(std::size_t n, const chainpars &par): mm2(0.0), pars(&par) { nmax=0; allocate(n);}

  chain(const chainpars &par): mm2(0.0), pars(&par), F_Pmod(false), F_Porigin(1), F_mm2s(false), num(0), nmax(0) {}

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
    F_mm2s=false;
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
    F_mm2s=false;
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
    memset(is_checked,false,num*sizeof(bool));
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
      for (std::size_t i=0;i<num;i++) {
        for (std::size_t j=i+1;j<num;j++) {
          mm2 += p[i].getMass() * p[j].getMass();
        }
      }
      mm2 /= num * (num - 1) / 2;
      F_mm2s = true;
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

      //Acceleration==========================================//
      memset(acc[lj],0,3*sizeof(double));
      
      //dw/dr ================================//
      memset(dWdr[lj],0,3*sizeof(double));
      
      for (std::size_t k=0;k<num;k++) {
        if(k==j) continue;
        std::size_t lk = list[k];
        const particle *pk= &p[lk];
        const double* xj = pj->getPos();
        double3 xjk;

        if(k==j+1) memcpy(xjk,X[j],3*sizeof(double));
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
          memset(At,0,3*sizeof(double));
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

  /* calc_dt
     function: calculate physical time step dt based on ds
     argument:  ds:  step size s (not physical time step) 
     return:    dt:  physical integration time step for X
  */
  double calc_dt(const double ds) {
    // determine the physical time step
    double dt = ds / (pars->alpha * (Ekin + B) + pars->beta * w + pars->gamma);
    if (std::abs(dt) < pars->dtmin) {
      std::cerr<<"Warning!: physical time step too small: "<<dt<<std::endl;
      abort();
    }
    return dt;
  }
  
  /* step_forward_X
     function: one step integration of X 
     argument: ds: step size s (not physical time step), it is modifed if final step reach
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
     argument: s: step size s (not physical time step)
     return: dt: time step for V (for B and w integration)
   */
  double step_forward_V(const double s) {
    // determine velocity integration time step
    double dt = s / (pars->alpha * Pot + pars->beta * W + pars->gamma);

    // step forward V
    for (std::size_t i=0;i<num-1;i++) {
      std::size_t k = list[i];
      std::size_t k1 = list[i+1];
      V[i][0] += dt * (acc[k1][0]-acc[k][0]);
      V[i][1] += dt * (acc[k1][1]-acc[k][1]);
      V[i][2] += dt * (acc[k1][2]-acc[k][2]);
    }

    return dt;
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
        memset(pf[i],0,3*sizeof(double));
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
              memcpy(xk,cj->p[k].getPos(),3*sizeof(double));
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
      memset(pf,0,3*num*sizeof(double));
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
    memcpy(roldlink,rlink,num*sizeof(std::size_t));

    // backup previous link
    std::size_t* listbk = new std::size_t[num];
    memcpy(listbk,list,num*sizeof(std::size_t));

    // backup current X
    double3* Xbk = new double3[num-1];
    memcpy(Xbk,X,(num-1)*3*sizeof(double));
    
    // backup current V
    double3* Vbk = new double3[num-1];
    memcpy(Vbk,V,(num-1)*3*sizeof(double));

    // create mask to avoid dup. check;
    bool* mask = new bool[num];
    memset(mask,false,num*sizeof(bool));

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

  /* Romberg_recursion_formula
     function: Using Romberg function: T_ik = T_i,k-1 + (T_i,k-1 - T_i-1,k-1) / [(h_i-k/h_i)^2 -1]
     argument: ti1k1: t_i-1,k-1
               tik1: t_i,k-1
               hr: h_i-k/h_i
     return:   T_ik
  */
  double Romberg_recursion_formula(const double ti1k1, const double tik1, const double hr) {
    if (hr==1) {
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
  double Rational_recursion_formula(const double ti1k2, const double ti1k1, const double tik1, const double hr) {
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
    pert_force();
    
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

  /* Leapfrog_step_forward
     function: integration n steps, particles' position and velocity will be updated in the center-of-mass frame
     argument: s: whole step size
               n: number of step division, real step size = (s/n), for Leapfrog integration, it is X(s/2n)V(s/n)X(s/n)V(s/n)..X(s/2n)
               toff: ending physical time (negative means no ending time and integration finishing after n steps)
               force: external force (not perturber forces which are calculated in pert_force)
               check_flag: 2: check link every step; 1: check link at then end; 0 :no check
               recur_flag: flag to determine whether to resolve sub-chain particles for force calculations. notice this require the sub-chain to be integrated to current physical time. Thus is this a recursion call (tree-recusion-integration)
///               upforce: void (const particle * p, const particle *pext, double3* force). function to calculate force based on p and pext, return to force
//// ext_force<particle> upforce) 
  */             
  void Leapfrog_step_forward(const double s, const int n, const double toff=-1.0, const double3* force=NULL, int check_flag=1, bool recur_flag=false) {
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
    if (pars->m_smooth&&!F_mm2s) {
      std::cerr<<"Error: smooth mass coefficients are used, but averaged mass coefficient mm2 is not calculated!\n";
      abort();
    }

    double ds = s/double(n);
    double3* ave_v=new double3[num];  // average velocity
    bool fpf = false; // perturber force indicator
    bool finalstep = false; 
    const int np = pext.getN();
    if (np>0) fpf = true;

    //initial step counter
    int i=1;

    //integration loop
    do {
      // half step forward for t (dependence: Ekin, B, w)
      double dt = calc_dt(ds*0.5);

      // if new time is larger than toff, cut it to fit toff.
      if (toff>0&&t+2*dt>toff) {
        dt = toff-t;
        ds = dt/(pars->alpha * (Ekin + B) + pars->beta * w + pars->gamma);
        dt = 0.5*dt; //get half dt
#ifdef DEBUG
        std::cerr<<"Reach time ending, current t="<<t<<"  toff="<<toff<<"  dt="<<dt<<"  ds="<<ds<<"  diff="<<toff-t<<std::endl;
#endif
      }
      // if toff not set, use n as finalstep criterion
      if (toff<0&&i==n) finalstep = true;

      // step_forward time
      t += dt;

      // recursive integration-----------------------------------------------//
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
      calc_rAPW(force, recur_flag); 

      // update chain list order if necessary, update list, X, V (dependence: p.x, X, V)
      if (num>2&&check_flag==2) update_link();

      // Step forward V and get time step dt(V) (dependence: V, Pot, A, W)
      double dvt = step_forward_V(ds);
      
      // Get averaged velocity, update p.x, p.v, ave_v (dependence: X, V)
      resolve_XV(ave_v);

      // forward B and w (dependence: dt(V), ave_v, p.m, p.v, force, dWdr, pf)
      step_forward_Bw(dvt,ave_v,force,fpf);

      // Calcuale Kinetic energy (dependence: p.m, p.v)
      calc_Ekin();

      // step forward for t (dependence: Ekin, B, w)
      dt = calc_dt(ds*0.5);
      t += dt;

      // check whether time match toff
      if (std::abs(toff - t)<pars->dterr) finalstep = true;
      
      // step forward for X (dependence: X, V)
      step_forward_X(dt);

      // reverst negative ds;
      if (ds<0) ds = -ds;
      
      i++;
    } while (!finalstep);

   // resolve X at last, update p.x (dependence: X)
    resolve_X();

    // Update rjk, A, Pot, dWdr, W (notice A will be incorrect since pf is not updated)
    calc_rAPW(force, recur_flag); 
    
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
   */
  double extrapolation_integration(const double ds, const double toff=-1.0, const double3* force=NULL) {
    // get parameters
    const double error = pars->exp_error;
    const std::size_t itermax = pars->exp_itermax;
    const int methods = pars->exp_methods;
    const std::size_t *step = pars->step;
    
#ifdef TIME_PROFILE
    profile.t_ep -= get_wtime();
#endif
    // array size indicator for relative position and velocity
    const std::size_t nrel = num-1;
    const std::size_t dsize = (2*nrel+1)*3;
    
    // for storage
    // for convenient, the data are storaged in one array with size (2*nrel+1)*3, which represents t, B, w, X[3][nrel], V[3][nrel]
    double d0[dsize],dtemp[dsize],dn[itermax][dsize];
    double Ekin0;

    // for error check
    double3 CX,CXN;
    double cxerr=error+1.0;
    double eerr=error+1.0;
    double cxerr0=cxerr+1.0;
    double eerr0=eerr+1.0;
    // error for step estimation
    double werrmax=std::numeric_limits<double>::max();

    // backup initial values
    d0[0] = t;
    d0[1] = B;
    d0[2] = w;
    memcpy(&d0[3],X,3*nrel*sizeof(double));
    memcpy(&d0[3+nrel*3],V,3*nrel*sizeof(double));
    Ekin0 = Ekin;

    // new step
    double dsn = 1.0;
    
    // for case of toff criterion
    /*double ds = s;
        if (toff>0) {
      // calculate ds from dt;
      ds = (toff-t)/(pars->alpha * (Ekin + B) + pars->beta * w + pars->gamma);
#ifdef DEBUG
      std::cerr<<"Extra-init, current t="<<t<<"  toff="<<toff<<"  dt="<<toff-t<<"  ds="<<ds<<std::endl;
#endif
} */     
    
    std::size_t intcount = 0; // iteration counter
    
    // first step
    Leapfrog_step_forward(ds,step[0],-1.0,force,0);

    std::size_t itercount = step[0]; // iteration efforts count

    // relative position vector between first and last particle for phase error check
    memset(CX,0,3*sizeof(double));
    for (size_t i=0;i<num-1;i++) {
      CX[0] += X[i][0];
      CX[1] += X[i][1];
      CX[2] += X[i][2];
    }

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
      if (intcount == 0) {
        // storage the results
        dn[0][0] = t;
        dn[0][1] = B;
        dn[0][2] = w;
        memcpy(&dn[0][3],X,3*nrel*sizeof(double));
        memcpy(&dn[0][3+nrel*3],V,3*nrel*sizeof(double));
      }
      intcount++;
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

      // reset the initial data
      t = d0[0];
      B = d0[1];
      w = d0[2];
      Ekin = Ekin0;
      memcpy(X,&d0[3],3*nrel*sizeof(double));
      memcpy(V,&d0[3+nrel*3],3*nrel*sizeof(double));
      // reset velocity to get correct w
      if (pars->beta>0) resolve_V();

      Leapfrog_step_forward(ds,step[intcount],-1.0,force,0);

      // increase iteration counter
      itercount += step[intcount];
      
      if (methods==1) {
        // Using Romberg method
        // first step extrapolate, storage the temp data in last index of array
        /*
          T_n-1,0 [0]
                     -> T_n,1 [n]
          T_n,0 [Chain]
        */
        // H_i-k/H_i
        double hr = (double)step[intcount]/(double)step[intcount-1];

        dn[intcount][0] = Romberg_recursion_formula(dn[0][0],t,hr);
        dn[intcount][1] = Romberg_recursion_formula(dn[0][1],B,hr);
        dn[intcount][2] = Romberg_recursion_formula(dn[0][2],w,hr);
        for (std::size_t i=0;i<nrel;i++) {
          for (std::size_t k=0;k<3;k++) {
            dn[intcount][3*(i+1)+k] = Romberg_recursion_formula(dn[0][3*(i+1)+k],X[i][k],hr);
            dn[intcount][3*(nrel+i+1)+k] = Romberg_recursion_formula(dn[0][3*(nrel+i+1)+k],V[i][k],hr);
          }
        }
   
        // if step number > 1, need iteration
        if (intcount>1) {
          // update zero order
          // [0] << T_n,0
          dn[0][0] = t;
          dn[0][1] = B;
          dn[0][2] = w;
          memcpy(&dn[0][3],X,3*nrel*sizeof(double));
          memcpy(&dn[0][3+nrel*3],V,3*nrel*sizeof(double));

          // iteration to get final result
          for (std::size_t j=1; j<intcount; j++) {
            //templately storage new results to ttemp
            /*
              T_n-1,j [j]
                       -> T_n,j+1 [temp]
              T_n,j [n]
             */
            hr = (double)step[intcount]/(double)step[intcount-j-1];
            for (std::size_t k=0; k<dsize; k++) dtemp[k] = Romberg_recursion_formula(dn[j][k],dn[intcount][k],hr);
   
            //update j order
            /*
              [j] << T_n,j [n]
             */
            memcpy(dn[j],dn[intcount],dsize*sizeof(double));
   
            //shift temp data to index = intcount.
            /*
              [n] << T_n,j+1 [temp]
             */
            memcpy(dn[intcount],dtemp,dsize*sizeof(double));
          }
        }
      }
      else {
        // Using Rational interpolation method
        // additional template storage
        double d1[dsize];
        /*double tt1, Bt1, wt1;
          double3* Xt1 = new double3[nrel];
          double3* Vt1 = new double3[nrel];*/
        
        /*
           T_n-1,0 [0]
        0               -> T_n,1 [n]
           T_n,0 [Chain]
        */
        double hr = (double)step[intcount]/(double)step[intcount-1];
        dn[intcount][0] = Rational_recursion_formula(0,dn[0][0],t,hr);
        dn[intcount][1] = Rational_recursion_formula(0,dn[0][1],B,hr);
        dn[intcount][2] = Rational_recursion_formula(0,dn[0][2],w,hr);
        for (std::size_t i=0;i<nrel;i++) {
          for (std::size_t k=0;k<3;k++) {
            dn[intcount][3*(i+1)+k] = Rational_recursion_formula(0,dn[0][3*(i+1)+k],X[i][k],hr);
            dn[intcount][3*(i+nrel+1)+k] = Rational_recursion_formula(0,dn[0][3*(i+nrel+1)+k],V[i][k],hr);
          }
        }

        // if step number > 1, need iteration
        if (intcount>1) {
          // iteration to get final result
          for (std::size_t j=1; j<intcount; j++) {
            //templately storage new results to ttemp
            /*
                           T_n-1,j [j]
             T_n-1,j-1 [j-1]          -> T_n,j+1 [temp]
                           T_n,j [n]
             */
            hr = (double)step[intcount]/(double)step[intcount-j-1];
            for (std::size_t k=0; k<dsize; k++) dtemp[k] = Rational_recursion_formula(dn[j-1][k],dn[j][k],dn[intcount][k],hr);

            if (j==1) {
              // update zero order
              // [0] << T_n,0
              dn[0][0] = t;
              dn[0][1] = B;
              dn[0][2] = w;
              memcpy(&dn[0][3],X,3*nrel*sizeof(double));
              memcpy(&dn[0][3+nrel*3],V,3*nrel*sizeof(double));
            }
            // update j-1 order
            // [j-1] << T_n,j-1 [t1]
            else memcpy(dn[j-1],d1,dsize*sizeof(double));

            //storage previous extrapolation data in template position t1
            // [t1] << T_n,j [n]
            memcpy(d1,dn[intcount],dsize*sizeof(double));
   
            //shift temp data to index = intcount.
            // [n] << T_n,j+1 [temp]
            memcpy(dn[intcount],dtemp,dsize*sizeof(double));
          }
        }
      }

      // get error estimation
      double ermax=0;
      for (std::size_t k=0; k<dsize; k++) {
        double dik = dn[intcount][k];
        double di1k= dn[intcount-1][k];
        ermax = std::max(ermax, 2*(dik-di1k)/std::sqrt(dik*dik+di1k*di1k));
      }
      double dsfactor = std::pow(ermax/error,1/double(2*intcount+3));
      double werrn = (double)itercount * dsfactor;
      if (ermax>0&&werrn<werrmax) {
        werrmax = werrn;
        dsn = 1.0 / dsfactor;
        //dsn = std::min(dsn,0.9); // not larger than 1.0
        //dsn = std::max(dsn,pars->dtmin); // not too small
#ifdef DEBUG
        std::cerr<<"ERR factor update: iterindex="<<intcount<<"; modify factor="<<1.0/dsfactor<<"; ermax="<<ermax<<std::endl;
#endif
      }
      
      // set final results back to chain array
      t = dn[intcount][0];
      B = dn[intcount][1];
      w = dn[intcount][2];
      memcpy(X,&dn[intcount][3],3*nrel*sizeof(double));
      memcpy(V,&dn[intcount][3+nrel*3],3*nrel*sizeof(double));
 
      // resolve particle
      resolve_XV();
      // recalculate the energy
      calc_Ekin();
      // force, potential and W
      calc_rAPW(force);

      // phase error calculation
      memset(CXN,0,3*sizeof(double));
      for (size_t i=0;i<nrel;i++) {
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
      memcpy(CX,CXN,3*sizeof(double));
#ifdef DEBUG
      std::cerr<<std::setprecision(6)<<"Iter.= "<<intcount<<" Dep.= "<<step[intcount]<<" P-err.= "<<cxerr;
      std::cerr<<" E-err.="<<(Ekin-Pot+B)/B<<" B ="<<std::setprecision(12)<<B<<std::endl;
#endif 
    }

    // update chain link order
    if(num>2) update_link();

#ifdef TIME_PROFILE
    profile.initstep(itermax+1);
    profile.stepcount[intcount+1]++;
    profile.t_ep += get_wtime();
#endif

    if (intcount+1<pars->opt_iter) dsn = std::pow(ds, (2*intcount+3)/(double)(2*pars->opt_iter+1)-1);
    return dsn;
  }

  /* getTime
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

