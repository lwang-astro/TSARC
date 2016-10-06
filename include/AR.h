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
#endif

//declaration

template <class particle> class chain;
template <class particle> class chainlist;

// external force calculation function pointer
//template <class particle>
//struct ext_force
//{
//  typedef void (*type)(const particle *, const particle *, double3*);
//};

// Chain class
template <class particle>
class chain{
  typedef double double3[3];
  double3 *X;  // relative position
  double3 *V;  // relative velocity
  std::size_t *list;   // chain index list
  double3 *acc; // acceleration
  double3 *pf;  // perturber force
  double3 *dmdr; // \partial Omega/ \partial rk

  //paramenters=======================================//
  double Ekin;  //kinetic energy
  double Pot;   //potential
  double W;     //time transformation function

  //Intgrt value======================================//
  double t;     //time
  double w;    //time transformation parameter
  double B;    //Binding energy (time momentum)

  std::size_t num;      //total number of chain particles
  std::size_t nmax;     //maximum number 

  //time step integration parameter
  double alpha; // logH cofficient
  double beta;  // TTL cofficient
  double gamma; // constant

  //  mass coefficients parameter
  double epi; // smooth parameter
  double m2;  // averaged mass
  bool m_smooth; // whether to use smooth mass coefficients 

  //monitor flags
  bool p_mod;     // indicate whether particle list is modified (true: modified)
  bool p_origin; // indicate whether particle is shifted back to original frame (true: original frame: false: center-of-mass frame)
  bool m2_s;     // indicate whether m2 is calculated.
  bool abg_s;    // indicate whether abg is set.

#ifdef TIME_PROFILE
  double t_apw;    // APW calculation
  double t_uplink; // update link
  double t_lf;     // leap-frog
  double t_ep;     // extrapolation
  double t_pext;   // perturber force
#endif

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
  
  //initialization
  //  template <class particle>
  chain(std::size_t n): alpha(1.0), beta(0.0), gamma(0.0), epi(0.001), m2(0.0), m_smooth(true), abg_s(true) { nmax=0; allocate(n);}

  chain(): alpha(1.0), beta(0.0), gamma(0.0), epi(0.001), m2(0.0), m_smooth(true), p_mod(false), p_origin(true), abg_s(true), m2_s(false), num(0), nmax(0) {}

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
    dmdr=new double3[n];
    p.init(n);
    p_mod=false;
    p_origin=true;
    m2_s=false;
  }

  // clear function
  void clear() {
    if (nmax>0) {
      delete[] X;
      delete[] V;
      delete[] list;
      delete[] acc;
      delete[] pf;
      delete[] dmdr;
      num = 0;
      nmax = 0;
    }
    p_mod=false;
    p_origin=true;
    m2_s=false;
  }

  //destruction
  ~chain() {
    if (nmax>0) {
      delete[] X;
      delete[] V;
      delete[] list;
      delete[] acc;
      delete[] pf;
      delete[] dmdr;
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


  /* calc_m2
     function: Get averaged mass coefficients for calculating Wjk
   */
  void calc_m2() {
      // calcualte m'^2
      for (std::size_t i=0;i<num;i++) {
        for (std::size_t j=i+1;j<num;j++) {
          m2 += p[i].getMass() * p[j].getMass();
        }
      }
      m2 /= num * (num - 1) / 2;
      m2_s = true;
  }    

  
  /* calc_wij
     function: Get mass coefficients from particle p.
               To avoid frequent shift of data, the i,j order follow particle p (not chain list)
     argument: option: 0): Wjk = m'^2 if m_i*m_j < epi*m'^2; 0 otherwise; 
                          m'^2 = sum (i<j) m_i*m_j/(N(N-1)/2)  (for case alpha=1, beta!=0)
                       1): Wjk = m_i*m_j      (for case alhpa=gamma=0, beta=1 or alpha=1, beta=gamma=0)
               epi: for option=0
  */  
  //  template <class particle>
  double calc_Wjk(const std::size_t i, const std::size_t j) {
    if(abg_s==1) {
      // Wjk = m2
      if (p[i].getMass()* p[j].getMass()<epi*m2) return m2;
      // Wjk = 0
      else return 0;
    }
    else {
      // Wjk = m_i*m_j
      return p[i].getMass() * p[j].getMass();
    }
  }

  /* calc_rAPW
     function: Get distance Matrix, acceleration, dm/dr, potential and transformation parameter
               based on particle mass (obtain from p), current X & V and Wjk
               (notice the acceleration array index k and the distance matrix index j,k order follow particle p to avoid additional shift when chain list change)
     argument: force: external force (acceleration) for each particle, (not perturber forces)
  */
  void calc_rAPW (const double3 *force=NULL) {
#ifdef TIME_PROFILE
    t_apw -= get_wtime();
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
      memset(dmdr[lj],0,3*sizeof(double));
      
      for (std::size_t k=0;k<num;k++) {
        if(k==j) continue;
        std::size_t lk = list[k];
        const particle *pk= &p[lk];
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
          const double* xj = pj->getPos();
          xjk[0] = xk[0] - xj[0];
          xjk[1] = xk[1] - xj[1];
          xjk[2] = xk[2] - xj[2];
        }

        double rjk = std::sqrt(xjk[0]*xjk[0]+xjk[1]*xjk[1]+xjk[2]*xjk[2]);
        double wjk = calc_Wjk(lj,lk);

        if (k>j) {
          //Potential energy==================================//
          Pot_c += pk->getMass() * pj->getMass() / rjk;
          //Transformation coefficient========================//
          W_c += wjk / rjk;
        }
        
        //Acceleration======================================//
        double rjk3 = rjk*rjk*rjk;
        double mor3 = pk->getMass() / rjk3;
        acc[lj][0] += mor3 * xjk[0];
        acc[lj][1] += mor3 * xjk[1];
        acc[lj][2] += mor3 * xjk[2];

        //d W / d r=========================================//
        mor3 = wjk / rjk3;
        dmdr[lj][0] += mor3 * xjk[0];
        dmdr[lj][1] += mor3 * xjk[1];
        dmdr[lj][2] += mor3 * xjk[2];
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
    t_apw += get_wtime();
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

  /* step_forward_X
     function: one step integration of X and physical time
     argument: ds: step size s (not physical time step), it is modifed if final step reach
               dtmin: minimum physical time step criterion
               toff: ending physical time (negative means no ending point)
     return:  final step flag (if toff is reached)
   */
  bool step_forward_X(double &ds, const double dtmin=5.4e-20, const double toff=-1) {
    bool finalstep = false;
    // determine the physical time step
    double weight = (alpha * (Ekin + B) + beta * w + gamma);
    double dt = ds / weight;
    if (dt<dtmin) {
      std::cerr<<"Warning!: physical time step too small: "<<dt<<std::endl;
      abort();
    }
    // if new time is larger than toff, cut it to fit toff.
    if (toff>0.0&&(t+2*dt>toff)) {
      dt = 0.5*(toff-t);
      ds = dt/ weight;
      finalstep=true;
    }

    // step forward physical time
    t += dt;

    // step forward relative X
    for (std::size_t i=0;i<num-1;i++) {
      X[i][0] += dt * V[i][0];
      X[i][1] += dt * V[i][1];
      X[i][2] += dt * V[i][2];
    }
    return finalstep;
  }

  /* step_forward_V
     function: one step integration of V
     argument: s: step size s (not physical time step)
     return: dt: time step for V (for B and w integration)
   */
  double step_forward_V(const double s) {
    // determine velocity integration time step
    double dt = s / (alpha * Pot + beta * W + gamma);

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
    if (beta>0||force!=NULL||fpf) {
      for (std::size_t i=0;i<num;i++) {
        if (force!=NULL) {
          dB -= p[i].getMass() * ( ave_v[i][0] * (pf[i][0] + force[i][0]) 
                            + ave_v[i][1] * (pf[i][1] + force[i][1]) 
                            + ave_v[i][2] * (pf[i][2] + force[i][2]));
        }
        else {
          dB -= p[i].getMass() * ( ave_v[i][0] * pf[i][0] 
                            + ave_v[i][1] * pf[i][1] 
                            + ave_v[i][2] * pf[i][2]);
        }          
        if (beta>0) {
          dw += ( ave_v[i][0] * dmdr[i][0]
                + ave_v[i][1] * dmdr[i][1]
                + ave_v[i][2] * dmdr[i][2]);
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
    argument: force: return pertuber force to each particles
    return:   true: pertubers exist. false: no perturbers
  */
  bool pert_force() {
#ifdef TIME_PROFILE
    t_pext -= get_wtime();
#endif
    const int np = pext.getN();
    if (np>0) {
      for (std::size_t i=0;i<num;i++) {
        memset(pf[i],0,3*sizeof(double));
        for (std::size_t j=0;j<np;j++) {
          const double *ri = p[i].getPos();
          const double *rpj = pext[j].getPos();
          const double mpj = pext[j].getMass();
          double dx = ri[0] - rpj[0];
          double dy = ri[1] - rpj[1];
          double dz = ri[2] - rpj[2];

          double rij2 = dx*dx + dy*dy + dz*dz;
          double rij3 = std::sqrt(rij2)*rij2;
        
          pf[i][0] -= mpj * dx / rij3;
          pf[i][1] -= mpj * dy / rij3;
          pf[i][2] -= mpj * dz / rij3;
        }
      }
      return true;
    }
    else {
      memset(pf,0,3*num*sizeof(double));
      return false;
    }
#ifdef TIME_PROFILE
    t_pext += get_wtime();
#endif
  }
     
  /* update_link
     function: update chain link based on current relative distance matrix and V
     return:   if link is modified, return true
   */
  bool update_link(){
#ifdef TIME_PROFILE
    t_uplink -= get_wtime();
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
    t_uplink += get_wtime();
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
    for (std::size_t i=0;i<num;i++) {
      const double *ri = p[i].getPos();
      const double *rc = cm.getPos();
      p[i].setPos(ri[0] + rc[0],
                  ri[1] + rc[1],
                  ri[2] + rc[2]);
    }
  }

  /* center_shift_X==================================
     function: shift positions and velocities of N (num) particles (p) based on chain center-of-mass
     p: particle list (read data and return shifted list);
  */
  void center_shift_X() {
    for (std::size_t i=0;i<num;i++) {
      const double *ri = p[i].getPos();
      const double *rc = cm.getPos();
      p[i].setPos(ri[0] - rc[0],
                  ri[1] - rc[1],
                  ri[2] - rc[2]);
    }
  }


public:
  /* add particle
     function: add particle into chainlist p
     argument: particle a
   */
  void addP(particle &a) {
    if (!p_origin) std::cerr<<"Warning!: particle list are in the center-of-mass frame, dangerous to add new particles!\n";
    p.add(a);
    p_mod=true;
  }

  void addP(chain<particle> &a) {
    if (!p_origin) std::cerr<<"Warning!: particle list are in the center-of-mass frame, dangerous to add new particles!\n";
    p.add(a);
    p_mod=true;
  }
  
  void addP(const std::size_t n, particle a[]) {
    if (!p_origin) std::cerr<<"Warning!: particle list are in the center-of-mass frame, dangerous to add new particles!\n";
    p.add(n,a);
    p_mod=true;
  }

  /* remove particle
     function: remove particle from p
     argument: i: particle index
               option: true: shift last particle to current position (defaulted);
                       false: shift all right particle to left by one
  */
  void removeP(const std::size_t i, bool option=true) { p.remove(i,option); p_mod=true; }


  /* initPext
     function: allocate array for Pext list
     argument: n: number of perturbers
   */
  void initPext(const std::size_t n) {
    if (pext.getN()) {
      std::cerr<"Error: Perturber list is already initialized!\n";
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
  
  /* is_p_modified
     function: reture true if particle list is modifed, in this case, the chain may need to be initialized again
     return: true/false
  */
  bool is_p_modified() const { return p_mod; }

  /* is_p_origin
     function: return true if particle list are in original frame
     return: true/false
  */
  bool is_p_origin() const { return p_origin; }
  
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
    if (p_origin) {
      center_shift_init();
      p_origin=false;
    }
    else {
      std::cerr<<"Error: particles are not in original frame!\n";
      abort();
    }

    //set member relative position and velocity
    calc_XV();
    
    // check whether time and mass coefficients are set, if not, use defaulted.
    if (!abg_s) setPars();

    // if smooth mass coefficients are used, calculate m2;
    if (m_smooth) calc_m2();

    //set relative distance matrix, acceleration, potential and transformation parameter
    calc_rAPW(force);

    //Initial intgrt value t
    t = time;

    //kinetic energy
    calc_Ekin();

    // initial time step parameter
    B = Pot - Ekin;
    w = W;

    // set p_mod to false
    p_mod = false;

#ifdef TIME_PROFILE
    reset_tp();
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
    std::cerr<<std::setw(WIDTH)<<Ekin<<std::setw(WIDTH)<<Pot<<std::setw(WIDTH)<<B+Ekin-Pot<<std::setw(WIDTH)<<B<<std::setw(WIDTH)<<w<<std::endl;
    #endif*/
    
  }


  /* center_shift_inverse==================================
     function: shift back to original positions and velocities of N (num) particles (p) based on their center-of-mass
     p: particle list (read data and return shifted list);
  */
  void center_shift_inverse() {
    if (p_origin) {
      std::cerr<<"Warning: particles are already in original frame!\n";
    }
    else {
      for (std::size_t i=0;i<num;i++) {
        const double *ri = p[i].getPos();
        const double *vi = p[i].getVel();
        const double *rc = cm.getPos();
        const double *vc = cm.getVel();
        p[i].setPos(ri[0] + rc[0],
                    ri[1] + rc[1],
                    ri[2] + rc[2]);
        p[i].setVel(vi[0] + vc[0],
                    vi[1] + vc[1],
                    vi[2] + vc[2]);
      }
      p_origin = true;
    }
  }

  /* center_shift==================================
     function: shift positions and velocities of N (num) particles (p) based on their center-of-mass
     p: particle list (read data and return shifted list);
  */
  void center_shift() {
    if (p_origin) {
      for (std::size_t i=0;i<num;i++) {
        const double *ri = p[i].getPos();
        const double *vi = p[i].getVel();
        const double *rc = cm.getPos();
        const double *vc = cm.getVel();
        p[i].setPos(ri[0] - rc[0],
                    ri[1] - rc[1],
                    ri[2] - rc[2]);
        p[i].setVel(vi[0] - vc[0],
                    vi[1] - vc[1],
                    vi[2] - vc[2]);
      }
      p_origin = false;
    }
    else {
      std::cerr<<"Warning: particles are already in the center-of-mass frame!\n";
    }      
  }

  /* setPars
     function: set time step integration parameter alpha, beta, gamma and mass coefficients parameters epi and smooth_flag
     argument: a: alpha (logH)
               b: beta  (TTL)
               g: gamma (constant)
               e: epi   (smooth parameter)
               option: m_smooth (whether use smooth mass coefficients)
               
  */
  void setPars(const double a=1.0, const double b=0.0, const double g=0.0, const double e=0.001, const bool mm=true) {
    alpha = a;
    beta = b;
    gamma = g;
    epi = e;
    m_smooth = mm;
    // safety check
    if (alpha==0&&beta==0&&gamma==0) {
      std::cerr<<"Error: alpha, beta and gamma cannot be all zero!\n";
      abort();
    }
    if (epi==0&&m_smooth) {
      std::cerr<<"Error: smooth mass coefficients are used, but smooth coefficient epi is zero!\n";
      abort();
    }
    abg_s = true;
  }


  /* Leapfrog_step_forward
     function: integration n steps, particles' position and velocity will be updated in the center-of-mass frame
     argument: s: whole step size
               n: number of step division, real step size = (s/n), for Leapfrog integration, it is X(s/2n)V(s/n)X(s/n)V(s/n)..X(s/2n)
               force: external force (not perturber forces which are calculated in pert_force)
               toff: ending physical time (negative means no ending time and integration finishing after n steps)
               dtmin: minimum physical time step criterion
               check_flag: 2: check link every step; 1: check link at then end; 0 :no check
///               force: external force array
///               upforce: void (const particle * p, const particle *pext, double3* force). function to calculate force based on p and pext, return to force
//// ext_force<particle> upforce) 
  */             
  void Leapfrog_step_forward(const double s, const int n, const double3* force=NULL, const double toff=-1.0, const double dtmin=5.4e-20, int check_flag=1) {
#ifdef TIME_PROFILE
    t_lf -= get_wtime();
#endif
    // Error check
    if (n<=0) {
      std::cerr<<"Error: step number shound be positive, current number is "<<n<<std::endl;
      abort();
    }
    if (s<0) {
      std::cerr<<"Error: step size should be positive, current value is "<<s<<std::endl;
      abort();
    }
    if (p_origin) {
      std::cerr<<"Error: particles are not in the center-of-mass frame, the integration can be dangerous!"<<std::endl;
      abort();
    }
    if (p_mod) {
      std::cerr<<"Error: particles are modified, initialization required!"<<std::endl;
      abort();
    }
    if (m_smooth&&!m2_s) {
      std::cerr<<"Error: smooth mass coefficients are used, but averaged mass coefficient m2 is not calculated!\n";
      abort();
    }

    double ds = s/double(n);
    double3* ave_v=new double3[num];  // average velocity
    bool fpf = false; // perturber force indicator
    const int np = pext.getN();
    if (np>0) fpf = true;

    // half step forward for X, t (dependence: Ekin, B, w, V)
    double hds = 0.5*ds;
    bool finalflag=step_forward_X(hds,dtmin,toff);

    // modify ds if finally step
    if (finalflag) ds = hds*2;

    //initial step counter
    int i=1;

    //integration loop
    do {
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
      
      // Update rjk, A, Pot, dmdr, W for half X (dependence: pf, force, p.getMass(), X, Wjk)
      calc_rAPW(force); 

      // update chain list order if necessary, update list, X, V (dependence: rjk, V)
      if (num>2&&check_flag==2) update_link();

      // Step forward V and get time step dt(V)
      double dvt = step_forward_V(ds);
      
      // Get averaged velocity, update p.x, p.v, ave_v (dependence: X, V)
      resolve_XV(ave_v);

      // forward B and w (dependence: dt(V), ave_v, p.getMass(), p.v, force, dmdr, pf)
      step_forward_Bw(dvt,ave_v,force,fpf);

      // Calcuale Kinetic energy (dependence: p.getMass(), p.v)
      calc_Ekin();

      // step forward for X (dependence: Ekin, B, w, V)
      if((toff<0&&i==n)||finalflag) {
        step_forward_X(hds,0.0);
        break;
      }
      else step_forward_X(ds,dtmin);
      
      i++;
    } while (true);

   // resolve X at last, update p.x (dependence: X)
    resolve_X();

    // Update rjk, A, Pot, dmdr, W (notice A will be incorrect since pf is not updated)
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
    t_lf += get_wtime();
#endif
  }


  /* extrapolation_integration
     function: extrapolation method to get accurate integration
     argument: s: whole step size
               error: relative error requirement for extrapolation
               itermax: maximum times for iteration.
               methods: 1: Romberg method; others: Rational interpolation method
               sequences: 1: even sequence {h, h/2, h/4, h/8 ...}; 2: Bulirsch & Stoer sequence {h, h/2, h/3, h/4, h/6, h/8 ...}
               force: external force (not perturber forces which are calculated in pert_force)
               dtmin: minimum physical time step
     return:   iteration count (start from 1, which means no iteration)
   */
  int extrapolation_integration(const double s, const double error=1E-8, const std::size_t itermax=10, const int methods=1, const int sequences=2, const double3* force=NULL, const double dtmin=5.4e-20) {
#ifdef TIME_PROFILE
    t_ep -= get_wtime();
#endif
    // array size indicator for relative position and velocity
    const std::size_t nrel = num-1;
    
    // for storage
    double   t0,ttemp,t1[itermax];
    double   B0,Btemp,B1[itermax];
    double   w0,wtemp,w1[itermax];
    double   Ekin0;

    double3* X1[itermax];
    double3* V1[itermax];
    for (std::size_t i=0;i<itermax;i++) {
      X1[i] = new double3[nrel];
      V1[i] = new double3[nrel];
    }
    double3* X0 = new double3[nrel];
    double3* V0 = new double3[nrel];
    double3* Xtemp = new double3[nrel];
    double3* Vtemp = new double3[nrel];

    // for error check
    double3 CX,CXN;
    //    double3 CXN;
    double cxerr=error+1.0;
    double eerr=error+1.0;
    double cxerr0=cxerr+1.0;
    double eerr0=eerr+1.0;

    // backup initial values
    t0 = t;
    B0 = B;
    w0 = w;
    Ekin0 = Ekin;
    memcpy(X0,X,3*nrel*sizeof(double));
    memcpy(V0,V,3*nrel*sizeof(double));
    //memcpy(p0,p,num*sizeof(particle));

    // first step
    Leapfrog_step_forward(s,1,force,-1.0,dtmin,0);

    // relative position vector between first and last particle for phase error check
    memset(CX,0,3*sizeof(double));
    for (size_t i=0;i<num-1;i++) {
      CX[0] += X[i][0];
      CX[1] += X[i][1];
      CX[2] += X[i][2];
    }

    std::size_t intcount = 0; // iteration counter
    std::size_t step[itermax]; //substep size
    step[0] = 1;
    // for B.S. 
    std::size_t stepeven = 2; 
    std::size_t stepodd = 3;
    
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
        t1[0] = t;
        B1[0] = B;
        w1[0] = w;
        memcpy(X1[0],X,3*nrel*sizeof(double));
        memcpy(V1[0],V,3*nrel*sizeof(double));
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

      if (sequences==1) {
        //even sequence
        step[intcount] = 2*step[intcount-1];
      }
      else {
        // B.S. sequence
        if (intcount%2) {
          step[intcount] = stepeven;
          stepeven = stepeven*2;
        }
        else {
          step[intcount] = stepodd;
          stepodd = stepodd*2;
        }
      }
      
      // reset the initial data
      t = t0;
      B = B0;
      w = w0;
      Ekin = Ekin0;
      memcpy(X,X0,3*nrel*sizeof(double));
      memcpy(V,V0,3*nrel*sizeof(double));
      // reset velocity to get correct w
      if (beta>0) resolve_V();

      Leapfrog_step_forward(s,step[intcount],force,-1.0,dtmin,0);

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
        
        t1[intcount] = Romberg_recursion_formula(t1[0],t,hr);
        B1[intcount] = Romberg_recursion_formula(B1[0],B,hr);
        w1[intcount] = Romberg_recursion_formula(w1[0],w,hr);
        for (std::size_t i=0;i<nrel;i++) {
          for (std::size_t k=0;k<3;k++) {
            X1[intcount][i][k] = Romberg_recursion_formula(X1[0][i][k],X[i][k],hr);
            V1[intcount][i][k] = Romberg_recursion_formula(V1[0][i][k],V[i][k],hr);
          }
        }
   
        // if step number > 1, need iteration
        if (intcount>1) {
          // update zero order
          // [0] << T_n,0
          t1[0] = t;
          B1[0] = B;
          w1[0] = w;
          memcpy(X1[0],X,3*nrel*sizeof(double));
          memcpy(V1[0],V,3*nrel*sizeof(double));

          // iteration to get final result
          for (std::size_t j=1; j<intcount; j++) {
            //templately storage new results to ttemp
            /*
              T_n-1,j [j]
                       -> T_n,j+1 [temp]
              T_n,j [n]
             */
            hr = (double)step[intcount]/(double)step[intcount-j-1];
            ttemp = Romberg_recursion_formula(t1[j],t1[intcount],hr);
            Btemp = Romberg_recursion_formula(B1[j],B1[intcount],hr);
            wtemp = Romberg_recursion_formula(w1[j],w1[intcount],hr);
            for (std::size_t i=0;i<nrel;i++) {
              for (std::size_t k=0;k<3;k++) {
                Xtemp[i][k] = Romberg_recursion_formula(X1[j][i][k],X1[intcount][i][k],hr);
                Vtemp[i][k] = Romberg_recursion_formula(V1[j][i][k],V1[intcount][i][k],hr);
              }
            }
   
            //update j order
            /*
              [j] << T_n,j [n]
             */
            t1[j] = t1[intcount];
            B1[j] = B1[intcount];
            w1[j] = w1[intcount];
            memcpy(X1[j],X1[intcount],3*nrel*sizeof(double));
            memcpy(V1[j],V1[intcount],3*nrel*sizeof(double));
   
            //shift temp data to index = intcount.
            /*
              [n] << T_n,j+1 [temp]
             */
            t1[intcount] = ttemp;
            B1[intcount] = Btemp;
            w1[intcount] = wtemp;
            memcpy(X1[intcount],Xtemp,3*nrel*sizeof(double));
            memcpy(V1[intcount],Vtemp,3*nrel*sizeof(double));
          }
        }
      }
      else {
        // Using Rational interpolation method
        // additional template storage
        double tt1, Bt1, wt1;
        double3* Xt1 = new double3[nrel];
        double3* Vt1 = new double3[nrel];
        
        /*
           T_n-1,0 [0]
        0               -> T_n,1 [n]
           T_n,0 [Chain]
        */
        double hr = (double)step[intcount]/(double)step[intcount-1];
        t1[intcount] = Rational_recursion_formula(0,t1[0],t,hr);
        B1[intcount] = Rational_recursion_formula(0,B1[0],B,hr);
        w1[intcount] = Rational_recursion_formula(0,w1[0],w,hr);
        for (std::size_t i=0;i<nrel;i++) {
          for (std::size_t k=0;k<3;k++) {
            X1[intcount][i][k] = Rational_recursion_formula(0,X1[0][i][k],X[i][k],hr);
            V1[intcount][i][k] = Rational_recursion_formula(0,V1[0][i][k],V[i][k],hr);
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
            ttemp = Rational_recursion_formula(t1[j-1],t1[j],t1[intcount],hr);
            Btemp = Rational_recursion_formula(B1[j-1],B1[j],B1[intcount],hr);
            wtemp = Rational_recursion_formula(w1[j-1],w1[j],w1[intcount],hr);
            for (std::size_t i=0;i<nrel;i++) {
              for (std::size_t k=0;k<3;k++) {
                Xtemp[i][k] = Rational_recursion_formula(X1[j-1][i][k],X1[j][i][k],X1[intcount][i][k],hr);
                Vtemp[i][k] = Rational_recursion_formula(V1[j-1][i][k],V1[j][i][k],V1[intcount][i][k],hr);
                //  if(Xtemp[i][k]!=Xtemp[i][k])  {
                //    std::cerr<<"j="<<j<<",i="<<i<<",k="<<k<<" Xik="<<Xtemp[i][k]<<" Xi1k2="<<X1[j-1][i][k]<<" Xi1k1="<<X1[j][i][k]<<" Xik1="<<X1[intcount][i][k]<<" hr="<<hr<<std::endl;
                //     abort();
                //  }
              }
            }

            if (j==1) {
              // update zero order
              // [0] << T_n,0
              t1[0] = t;
              B1[0] = B;
              w1[0] = w;
              memcpy(X1[0],X,3*nrel*sizeof(double));
              memcpy(V1[0],V,3*nrel*sizeof(double));
            }
            else {
              // update j-1 order
              // [j-1] << T_n,j-1 [t1]
              t1[j-1] = tt1;
              B1[j-1] = Bt1;
              w1[j-1] = wt1;
              memcpy(X1[j-1],Xt1,3*nrel*sizeof(double));
              memcpy(V1[j-1],Vt1,3*nrel*sizeof(double));
            }

            //storage previous extrapolation data in template position t1
            // [t1] << T_n,j [n]
            tt1 = t1[intcount];
            Bt1 = B1[intcount];
            wt1 = w1[intcount];
            memcpy(Xt1,X1[intcount],3*nrel*sizeof(double));
            memcpy(Vt1,V1[intcount],3*nrel*sizeof(double));
   
            //shift temp data to index = intcount.
            // [n] << T_n,j+1 [temp]
            t1[intcount] = ttemp;
            B1[intcount] = Btemp;
            w1[intcount] = wtemp;
            memcpy(X1[intcount],Xtemp,3*nrel*sizeof(double));
            memcpy(V1[intcount],Vtemp,3*nrel*sizeof(double));
          }
        }

        delete[] Xt1;
        delete[] Vt1;
      }

      // set final results back to chain array
      t = t1[intcount];
      B = B1[intcount];
      w = w1[intcount];
      memcpy(X,X1[intcount],3*nrel*sizeof(double));
      memcpy(V,V1[intcount],3*nrel*sizeof(double));
 
      // resolve particle
      resolve_XV();
      // recalculate the energy
      calc_Ekin();
      // force, potential and W
      calc_rAPW();

      // phase error calculation
      memset(CXN,0,3*sizeof(double));
      for (size_t i=0;i<num-1;i++) {
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
      std::cerr<<std::setprecision(18)<<"Iteration: "<<intcount<<" Step division = "<<step[intcount]<<" error = "<<cxerr;
      std::cerr<<" energy error "<<Ekin-Pot+B<<std::endl;
#endif 
    }

    delete[] Xtemp;
    delete[] Vtemp;
    delete[] X0;
    delete[] V0;
    for (std::size_t i=0;i<itermax;i++) {
      delete[] X1[i];
      delete[] V1[i];
    }
    
    // update chain link order
    if(num>2) update_link();

#ifdef TIME_PROFILE
    t_ep += get_wtime();
#endif
    
    return intcount+1;
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
  
#ifdef TIME_PROFILE
  void reset_tp(){
    t_apw=0.0;
    t_uplink=0.0;
    t_lf=0.0;
    t_ep=0.0;
    t_pext=0.0;
  }

  // calc_APW
  double getTP_apw() const{
    return t_apw;
  }

  //update link
  double getTP_uplink() const{
    return t_uplink;
  }

  // leap-frog
  double getTP_lf() const{
    return t_lf;
  }

  // extrapolation
  double getTP_ep() const{
    return t_ep;
  }

  // perturber force
  double getTP_pext() const{
    return t_pext;
  }
#endif


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
      for (std::size_t i=0;i<num;i++) std::cerr<<std::setw(width)<<dmdr[i][k];
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
             <<"alpha: "<<std::setw(width)<<alpha<<std::endl
             <<"beta: "<<std::setw(width)<<beta<<std::endl
             <<"gamma: "<<std::setw(width)<<gamma<<std::endl;
  }
#endif  
};

// Generalized list to store chain and Particle members
template <class particle>
class chainlist{
  std::size_t num;
  std::size_t nmax;
  bool* cflag; // flag to indicater whether it is chain (true: chain; false: Particle)
  void** p;

public:
  // initialization
  chainlist(): num(0), nmax(0) {};
  
  chainlist(const std::size_t n) { init(n); }

  void init(const std::size_t n) {
    num = 0;
    nmax = n;
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
  void remove(std::size_t i, bool option=true) {
    if (option) {
      if (i<num-1) {
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
  particle &operator [](std::size_t i){
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

};

