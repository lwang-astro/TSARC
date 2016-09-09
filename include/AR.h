#pragma once

#include <cassert>
#include <iostream>
#include <string.h>
#include <cmath>
#include <stdlib.h>
#include <iomanip>

#ifndef NAN_CHECK
#define NAN_CHECK(val) assert((val) == (val));
#endif

#ifdef DEBUG
#define WIDTH 10
#endif

typedef double double3[3];
typedef double double4[4];

// external force calculation function pointer
//template <class particle>
//struct ext_force
//{
//  typedef void (*type)(const particle *, const particle *, double3*);
//};

//basic particle structure==========================//
class Particle{
public:
  double3 pos, vel;
  double mass;

  //instruction=======================================//
  Particle() {} 
  Particle(double m, double r[3], double v[3]) {
    load(m,r[0],r[1],r[2],v[0],v[1],v[2]);
  }
  Particle(double m, double rx, double ry, double rz, double vx, double vy, double vz) {
    load(m,rx,ry,rz,vx,vy,vz);
  }

  // Get position
  double3* getPos() {
    return &pos;
  }

  // Get velocity
  double3* getVel() {
    return &vel;
  }
  
  //load data=========================================//
  void load(double m, double rx, double ry, double rz, double vx, double vy, double vz){
    NAN_CHECK(m);
    NAN_CHECK(rx);
    NAN_CHECK(ry);
    NAN_CHECK(rz);
    NAN_CHECK(vx);
    NAN_CHECK(vy);
    NAN_CHECK(vz);

    mass=m;
    pos[0]=rx;
    pos[1]=ry;
    pos[2]=rz;
    vel[0]=vx;
    vel[1]=vy;
    vel[2]=vz;
  }

};

// Chain class
template <class particle>
class chain{
private:
  double3 *X;  // relative position
  double3 *V;  // relative velocity
  int *list;   // chain index list
  double3 *acc; // acceleration
  double3 *pf;  // perturber force
  double3 *dmdr; // \partial Omega/ \partial rk
  double **Wjk; // Omega ij, mass coefficients
  double4 **rjk; // relative position matrix

  //paramenters=======================================//
  double Ekin;  //kinetic energy
  double Pot;   //potential
  double W;     //time transformation function

  //Intgrt value======================================//
  double t;     //time
  double w;    //time transformation parameter
  double B;    //Binding energy (time momentum)

  int num;      //total number of chain particles

  //time step integration parameter
  double alpha; 
  double beta;
  double gamma;
  double dtmin; // minimum time step for problem checking.

public:

  //center mass=======================================//
  particle cm;
  //  double cm.mass;  //total mass
  //  double3 cm.pos;   // position center
  //  double3 cm.vel;   // velocity center

  //initialization
  //  template <class particle>
  chain(int n): num(n) {
    X=new double3[n-1];
    V=new double3[n-1];
    list=new int[n];
    acc=new double3[n];
    pf=new double3[n];
    dmdr=new double3[n];
    Wjk=new double*[n];
    rjk=new double4*[n];
    for (int i=0;i<n;i++) {
      Wjk[i]=new double[n];
      rjk[i]=new double4[n];
    }
    dtmin=0;
  }

private:
  /* update_num
     function: to update the num due to current particle list number
     argument: n: particle number
  */
  void update_num(const int n) {
    if (n>num) {
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
     argument: p: particle list arry
  */
  //  template <class particle>
  void generate_list(const particle *p) {
    bool *is_checked=new bool[num];
    memset(is_checked,false,num*sizeof(bool));
    int inext=0;
    for (int i=0; i<num; i++) {
      //initial rjk=======================================//
      //mark checked particle=============================//
      is_checked[inext] = true;

      //initial chain_mem=================================//
      list[i]=inext;
      int inow=inext;
    
      //make chain========================================//
      double rmin;
      bool first=true;
      for (int j=1; j<num; j++) {
        if(is_checked[j]) continue;
        double dx=p[j].pos[0] - p[inow].pos[0];
        double dy=p[j].pos[1] - p[inow].pos[1];
        double dz=p[j].pos[2] - p[inow].pos[2];
        double dr2=dx*dx + dy*dy + dz*dz;
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

  /* set_center
     function: set center-of-mass information
     argument: p: center-of-mass particle 
   */
  //  template <class particle>
  void set_center(const particle &p) {
    memcpy(cm.pos,p.pos,3*sizeof(double));
    memcpy(cm.vel,p.vel,3*sizeof(double));
    cm.mass = p.mass;
  }
  
  /* calc_XV
     function: get chain member relative position X and velocity V based on current list
     argument: p: particle list
  */
  //  template <class particle>
  void calc_XV(const particle *p) {
    for (int i=0;i<num-1;i++) {
      X[i][0]=p[list[i+1]].pos[0] - p[list[i]].pos[0];
      X[i][1]=p[list[i+1]].pos[1] - p[list[i]].pos[1];
      X[i][2]=p[list[i+1]].pos[2] - p[list[i]].pos[2];

      V[i][0]=p[list[i+1]].vel[0] - p[list[i]].vel[0];
      V[i][1]=p[list[i+1]].vel[1] - p[list[i]].vel[1];
      V[i][2]=p[list[i+1]].vel[2] - p[list[i]].vel[2];
    }
  }

  /* calc_wij
     function: Get mass coefficients from particle p.
               To avoid frequent shift of data, the i,j order follow particle p (not chain list)
     argument: p: particle list
     option: 0): Wjk = m'^2 if m_i*m_j < epi*m'^2; 0 otherwise; 
     m'^2 = sum (i<j) m_i*m_j/(N(N-1)/2)  (for case alpha=1, beta!=0)
     1): Wjk = m_i*m_j      (for case alhpa=gamma=0, beta=1 or alpha=1, beta=gamma=0)
     epi: for option=0
  */  
  //  template <class particle>
  void calc_Wjk(const particle *p, const int option=0, const double epi=0.001) {
    if(option==0) {
      // calcualte m'^2
      double m2=0.;
      for (int i=0;i<num;i++) {
        for (int j=i+1;j<num;j++) {
          m2 += p[i].mass * p[j].mass;
        }
      }
      m2 /= num * (num - 1) / 2;

      // set Wjk
      for (int i=0;i<num;i++) {
        for (int j=0;j<num;j++) {
          if (p[i].mass*p[j].mass<epi*m2) Wjk[i][j] = m2;
          else Wjk[i][j]=0;
        }
      }
    }
    else {
      // Wjk = m_i*m_j
      for (int i=0;i<num;i++) {
        for (int j=0;j<num;j++) {
          Wjk[i][j] = p[i].mass * p[j].mass;
        }
      }
    }
  }

  /* calc_rAPW
     function: Get distance Matrix, acceleration, dm/dr, potential and transformation parameter
               based on particle mass (obtain from p), current X & V and Wjk
               (notice the acceleration array index k and the distance matrix index j,k order follow particle p to avoid additional shift when chain list change)
     argument: p: particle list (only mass needed)
               force: external force (acceleration) for each particle, (not perturber forces)
  */
  //  template <class particle>
  void calc_rAPW (const particle *p, const double3 *force=NULL) {
    // template xjk vector
    double3 xjk={0};
    // reset potential and transformation parameter
    Pot  = 0.0;
    W  = 0.0;
    // Loop all particles 
    for (int j=0;j<num;j++) {
      int lj = list[j];
      const particle *pj= &p[lj];

      //Acceleration==========================================//
      memset(acc[lj],0,3*sizeof(double));
      
      //dw/dr ================================//
      memset(dmdr[lj],0,3*sizeof(double));
      
      for (int k=0;k<num;k++) {
        int lk = list[k];
        const particle *pk= &p[lk];
        double *rljk= rjk[lj][lk];

        if(k==j) memset(rljk,0,4*sizeof(double));
        else if (k>j) {
          if(k==j+1) memcpy(xjk,X[j],3*sizeof(double)); // initial xjk;
          else { // update xjk
            xjk[0] += X[k-1][0];
            xjk[1] += X[k-1][1];
            xjk[2] += X[k-1][2];
          }
          memcpy(rljk,xjk,3*sizeof(double));
          // Suppressed old method
          // else if(k==j+2) {
          //   rljk[0] = X[j][0] + X[j+1][0];
          //   rljk[1] = X[j][1] + X[j+1][1];
          //   rljk[2] = X[j][2] + X[j+1][2];
          // }
          // else {
          //   rljk[0] = pk->x[0] - pj->x[0];
          //   rljk[1] = pk->x[1] - pj->x[1];
          //   rljk[2] = pk->x[2] - pj->x[2];

          rljk[3] = std::sqrt(rljk[0]*rljk[0]+rljk[1]*rljk[1]+rljk[2]*rljk[2]);
          //Potential energy==================================//
          Pot += pk->mass * pj->mass / rljk[3];
          //Transformation coefficient========================//
          W += Wjk[lj][lk] / rljk[3];
        }
        else {
          rljk[0] = -rjk[lk][lj][0];
          rljk[1] = -rjk[lk][lj][1];
          rljk[2] = -rjk[lk][lj][2];
          rljk[3] =  rjk[lk][lj][3];
        }
        
        if (k!=j) {
          //Acceleration======================================//
          double rljk3 = rljk[3]*rljk[3]*rljk[3];
          double mor3 = pk->mass / rljk3;
          acc[lj][0] += mor3 * rljk[0];
          acc[lj][1] += mor3 * rljk[1];
          acc[lj][2] += mor3 * rljk[2];

          //d W / d r=========================================//
          mor3 = Wjk[lj][lk] / rljk3;
          dmdr[lj][0] += mor3 * rljk[0];
          dmdr[lj][1] += mor3 * rljk[1];
          dmdr[lj][2] += mor3 * rljk[2];
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
  }

  /* calc_Ekin
     function: calculate kinetic energy
  */
  //  template <class particle>
  void calc_Ekin(const particle* p){
    Ekin = 0.0;
    for (int i=0; i<num; i++)
      Ekin += 0.5 * p[i].mass * (p[i].vel[0]*p[i].vel[0]+p[i].vel[1]*p[i].vel[1]+p[i].vel[2]*p[i].vel[2]);
  }

  /* calc_Pot
     function: calculate potential energy
  */
  //  template <class particle>
  void calc_Pot(const particle* p){
    Pot = 0.0;
    for (int i=0; i<num; i++) 
      for (int j=0; j<i; j++) 
        Pot += p[i].mass * p[j].mass / rjk[i][j][3];
  }
  

  /* step_forward_X
     function: one step integration of X and physical time
     argument: step size s (not physical time step)
   */
  void step_forward_X(const double s) {
    // determine the physical time step
    double dt = s / (alpha * (Ekin + B) + beta * w + gamma);
    if (dt<dtmin) {
      std::cerr<<"Warning!: physical time step too small: "<<dt<<std::endl;
      abort();
    }

    // step forward physical time
    t += dt;

    // step forward relative X
    for (int i=0;i<num-1;i++) {
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
    double dt = s / (alpha * Pot + beta * W + gamma);

    // step forward V
    for (int i=0;i<num-1;i++) {
      int k = list[i];
      int k1 = list[i+1];
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
  void step_forward_Bw(const double dt, const double3* ave_v, const particle* p, const double3* force, const bool fpf) {
    double dB = 0.0;
    double dw = 0.0;
    if (beta>0||force!=NULL||fpf) {
      for (int i=0;i<num;i++) {
        if (force!=NULL) {
          dB -= p[i].mass * ( ave_v[i][0] * (pf[i][0] + force[i][0]) 
                            + ave_v[i][1] * (pf[i][1] + force[i][1]) 
                            + ave_v[i][2] * (pf[i][2] + force[i][2]));
        }
        else {
          dB -= p[i].mass * ( ave_v[i][0] * pf[i][0] 
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
               The total mass of particles should be consistent with cm.mass. Otherwise update Chain.cm first.
     argument: ave_v: averaged velocity array (return values)
               p: particle lists (return new x, v with center-of-mass correction)
   */
  //  template <class particle>
  void resolve_XV(particle* p, double3* ave_v=NULL) {
    // backup old v
    if (ave_v!=NULL) {
      for (int i=0;i<num;i++) {
        ave_v[i][0] = p[i].vel[0];
        ave_v[i][1] = p[i].vel[1];
        ave_v[i][2] = p[i].vel[2];
      }
    }
    // resolve current V
    double3 vc={0};
    double3 xc={0};
    for (int i=0;i<num-1;i++) {
      int lk = list[i];
      int lkn = list[i+1];
      p[lkn].pos[0] = p[lk].pos[0] + X[i][0];
      p[lkn].pos[1] = p[lk].pos[1] + X[i][1];
      p[lkn].pos[2] = p[lk].pos[2] + X[i][2];
      p[lkn].vel[0] = p[lk].vel[0] + V[i][0];
      p[lkn].vel[1] = p[lk].vel[1] + V[i][1];
      p[lkn].vel[2] = p[lk].vel[2] + V[i][2];
      //center-of-mass position and velocity
      xc[0] += p[lkn].mass * p[lkn].pos[0];
      xc[1] += p[lkn].mass * p[lkn].pos[1];
      xc[2] += p[lkn].mass * p[lkn].pos[2];
      vc[0] += p[lkn].mass * p[lkn].vel[0];
      vc[1] += p[lkn].mass * p[lkn].vel[1];
      vc[2] += p[lkn].mass * p[lkn].vel[2];
      if (i==0) {
        xc[0] += p[lk].mass * p[lk].pos[0];
        xc[1] += p[lk].mass * p[lk].pos[1];
        xc[2] += p[lk].mass * p[lk].pos[2];
        vc[0] += p[lk].mass * p[lk].vel[0];
        vc[1] += p[lk].mass * p[lk].vel[1];
        vc[2] += p[lk].mass * p[lk].vel[2];
      }        
    }

    // calcualte center-of-mass position and velocity shift
    xc[0] /= cm.mass;
    xc[1] /= cm.mass;
    xc[2] /= cm.mass;
    vc[0] /= cm.mass;
    vc[1] /= cm.mass;
    vc[2] /= cm.mass;
    
    for (int i=0;i<num;i++) {
      // center-of-mass correction
      p[i].pos[0] -= xc[0];
      p[i].pos[1] -= xc[1];
      p[i].pos[2] -= xc[2];
      p[i].vel[0] -= vc[0];
      p[i].vel[1] -= vc[1];
      p[i].vel[2] -= vc[2];

      // calculate averaged velocities
      if (ave_v!=NULL) {
        ave_v[i][0] = 0.5 * (ave_v[i][0] + p[i].vel[0]);
        ave_v[i][1] = 0.5 * (ave_v[i][1] + p[i].vel[1]);
        ave_v[i][2] = 0.5 * (ave_v[i][2] + p[i].vel[2]);
      }
    }
  }

  /* resolve_X
     function: resolve relative X to physical x (center-of-mass frame)
               Notice the center-of-mass particle mass in Chain.cm is used.
               The total mass of particles should be consistent with cm.mass. Otherwise update Chain.cm first.
     argument: p: particle lists (return new x with center-of-mass correction)
   */
  //  template <class particle>
  void resolve_X(particle* p) {
    // resolve current X
    double3 xc={0};
    for (int i=0;i<num-1;i++) {
      int lk = list[i];
      int lkn = list[i+1];
      p[lkn].pos[0] = p[lk].pos[0] + X[i][0];
      p[lkn].pos[1] = p[lk].pos[1] + X[i][1];
      p[lkn].pos[2] = p[lk].pos[2] + X[i][2];
      //center-of-mass position and velocity
      xc[0] += p[lkn].mass * p[lkn].pos[0];
      xc[1] += p[lkn].mass * p[lkn].pos[1];
      xc[2] += p[lkn].mass * p[lkn].pos[2];
      if (i==0) {
        xc[0] += p[lk].mass * p[lk].pos[0];
        xc[1] += p[lk].mass * p[lk].pos[1];
        xc[2] += p[lk].mass * p[lk].pos[2];
      }
    }

    // calcualte center-of-mass position and velocity shift
    xc[0] /= cm.mass;
    xc[1] /= cm.mass;
    xc[2] /= cm.mass;
    
    for (int i=0;i<num;i++) {
      // center-of-mass correction
      p[i].pos[0] -= xc[0];
      p[i].pos[1] -= xc[1];
      p[i].pos[2] -= xc[2];
    }
  }

  /* pert_force
    function: get perturber force
    argument: p: particle list
              pext: pertuber particle list
              np: perturber number
              force: return pertuber force to each particles
    return:   true: pertubers exist. false: no perturbers
  */
  bool pert_force(const particle *p, const int np, const particle *pext) {
    if (np>0) {
      for (int i=0;i<num;i++) {
        memset(pf[i],0,3*sizeof(double));
        for (int j=0;j<np;j++) {
          double dx = p[i].pos[0] - pext[j].pos[0];
          double dy = p[i].pos[1] - pext[j].pos[1];
          double dz = p[i].pos[2] - pext[j].pos[2];

          double rij2 = dx*dx + dy*dy + dz*dz;
          double rij3 = std::sqrt(rij2)*rij2;
        
          pf[i][0] -= pext[j].mass * dx / rij3;
          pf[i][1] -= pext[j].mass * dy / rij3;
          pf[i][2] -= pext[j].mass * dz / rij3;
        }
      }
      return true;
    }
    else {
      memset(pf,0,3*num*sizeof(double));
      return false;
    }
  }
     
public:
  /* init
     function: initialization chain from particle p
     argument: time: current time
     n:    number of particles
     p:    particle list (particle will be shifted based on center-of-mass)
     pext: perturber particle list
     np:   perturber particle number
     force: external force
  */
  
  //  template <class particle>
  void init(const double time, const int n, particle *p, const int np=0, const particle *pext=NULL, const double3* force=NULL) {
    // adjust num to current n
    update_num(n);

    //Generate chain link list
    generate_list(p);

    // calculate perturber forces
    pert_force(p,np,pext);
    
    //set center-of-mass
    center_shift_init(p);

    //set member relative position and velocity
    calc_XV(p);
    
    //set mass coefficients Wjk
    if (gamma==0&&((alpha==1&&beta==0)||(alpha==0&&beta==1))) calc_Wjk(p,1);
    else calc_Wjk(p);

    //set relative distance matrix, acceleration, potential and transformation parameter
    calc_rAPW(p,force);

    //Initial intgrt value t
    t = time;

    //kinetic energy
    calc_Ekin(p);

    //Potential energy
    //calc_Pot(p);
    
    // initial time step parameter
    B = Pot - Ekin;
    w = B;
    /*#ifdef DEBUG
//    for (int i=0;i<num;i++) {
//      std::cout<<i<<" m"<<std::setw(WIDTH)<<p[i].mass<<std::setw(WIDTH)<<"x";
//      for (int k=0;k<3;k++) std::cout<<std::setw(WIDTH)<<p[i].pos[k];
//      std::cout<<std::setw(WIDTH)<<"v ";
//      for (int k=0;k<3;k++) std::cout<<std::setw(WIDTH)<<p[i].vel[k];
//      std::cout<<std::endl;
//    }
//    std::cout<<std::setw(WIDTH)<<"kin:"<<std::setw(WIDTH)<<0.5*(p[0].mass*p[1].mass/(p[0].mass+p[1].mass))*(V[0][0]*V[0][0]+V[0][1]*V[0][1]+V[0][2]*V[0][2])<<std::endl;
    std::cerr<<std::setw(WIDTH)<<t;
    for (int i=0;i<num-1;i++) {
      for (int k=0;k<3;k++) std::cerr<<std::setw(WIDTH)<<X[i][k];
      for (int k=0;k<3;k++) std::cerr<<std::setw(WIDTH)<<V[i][k];
    }
    std::cerr<<std::setw(WIDTH)<<Ekin<<std::setw(WIDTH)<<Pot<<std::setw(WIDTH)<<B+Ekin-Pot<<std::setw(WIDTH)<<B<<std::setw(WIDTH)<<w<<std::endl;
    #endif*/
  }

  /* update_link
     function: update chain link based on current relative distance matrix and V
   */
  void update_link(){
    // create reverse index of link
    int* rlink = new int[num];
    for (int i=0;i<num;i++) rlink[list[i]] = i;

    // backup previous link
    int* listbk = new int[num];
    memcpy(listbk,list,num*sizeof(double));
    
    // backup current V
    double3* Vbk = new double3[num];
    memcpy(Vbk,V,num*3*sizeof(double));

    // create mask to avoid dup. check;
    bool* mask = new bool[num];
    memset(mask,false,num*sizeof(bool));
    for (int k=0;k<num-1;k++) {
      int lk  = list[k];
      mask[lk] = true;
      int lkn = list[k+1];
      // possible new index
      int lku = lkn;
      double rmin = rjk[lk][lkn][3];
      for (int j=0;j<num;j++) {
        if (mask[j]||j==k) continue;
        if (rjk[lk][j][3]<rmin) {
          lku=j;
          rmin = rjk[lk][j][3];
        }
      }
      if (lku!=lkn) {
        // shift two index in the list
        list[rlink[lku]] = lkn;
        list[k+1] = lku;
        mask[lku] = true;
      }

      if (lk!=listbk[k]||lku!=lkn) {
        // update X from rjk
        memcpy(X[k], rjk[lk][lku], 3*sizeof(double));
        // update V
        // left boundary
        int rlk = rlink[lk];
        if (rlk<k) {
          for (int j=rlk;j<k;j++) {
            V[k][0] += Vbk[j][0];
            V[k][1] += Vbk[j][1];
            V[k][2] += Vbk[j][2];
          }
        }
        else if (rlk>k) {
          for (int j=k;j<rlk;j++) {
            V[k][0] -= Vbk[j][0];
            V[k][1] -= Vbk[j][1];
            V[k][2] -= Vbk[j][2];
          }
        }
        // right boundary
        rlk = rlink[lku];
        if (rlk<k+1) {
          for (int j=rlk;j<k+1;j++) {
            V[k][0] -= Vbk[j][0];
            V[k][1] -= Vbk[j][1];
            V[k][2] -= Vbk[j][2];
          }
        }
        else if (rlk>k+1) {
          for (int j=k+1;j<rlk;j++) {
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
    delete[] listbk;
    delete[] Vbk;
  }

  /* center_shift ==================================
     function: shift positions and velocities of N (num) particles (p) based on their center-of-mass, write center-of-mass particle to chain
     p: particle list (read data and return shifted list);
  */
  void center_shift_init(particle *p) {
    //center mass=======================================//
    memset(cm.pos,0,3*sizeof(double));
    memset(cm.vel,0,3*sizeof(double));
    cm.mass = 0.0;
    for (int i=0;i<num;i++) {
      cm.pos[0] += p[i].pos[0] * p[i].mass;
      cm.pos[1] += p[i].pos[1] * p[i].mass;
      cm.pos[2] += p[i].pos[2] * p[i].mass;

      cm.vel[0] += p[i].vel[0] * p[i].mass;
      cm.vel[1] += p[i].vel[1] * p[i].mass;
      cm.vel[2] += p[i].vel[2] * p[i].mass;

      cm.mass += p[i].mass;
    }
    cm.pos[0] /= cm.mass; 
    cm.pos[1] /= cm.mass; 
    cm.pos[2] /= cm.mass; 

    cm.vel[0] /= cm.mass; 
    cm.vel[1] /= cm.mass; 
    cm.vel[2] /= cm.mass; 

    // shifting
    for (int i=0;i<num;i++) {
      p[i].pos[0] -= cm.pos[0];
      p[i].pos[1] -= cm.pos[1];
      p[i].pos[2] -= cm.pos[2];
      p[i].vel[0] -= cm.vel[0];
      p[i].vel[1] -= cm.vel[1];
      p[i].vel[2] -= cm.vel[2];
    }
  }

  /* center_shift_inverse==================================
     function: shift back positions and velocities of N (num) particles (p) based on their center-of-mass, write center-of-mass particle to chain
     p: particle list (read data and return shifted list);
  */
  void center_shift_inverse(particle *p) {
    for (int i=0;i<num;i++) {
      p[i].pos[0] += cm.pos[0];
      p[i].pos[1] += cm.pos[1];
      p[i].pos[2] += cm.pos[2];
      p[i].vel[0] += cm.vel[0];
      p[i].vel[1] += cm.vel[1];
      p[i].vel[2] += cm.vel[2];
    }
  }

  /* center_shift==================================
     function: shift positions and velocities of N (num) particles (p) based on their center-of-mass, write center-of-mass particle to chain
     p: particle list (read data and return shifted list);
  */
  void center_shift(particle *p) {
    for (int i=0;i<num;i++) {
      p[i].pos[0] -= cm.pos[0];
      p[i].pos[1] -= cm.pos[1];
      p[i].pos[2] -= cm.pos[2];
      p[i].vel[0] -= cm.vel[0];
      p[i].vel[1] -= cm.vel[1];
      p[i].vel[2] -= cm.vel[2];
    }
  }

  /* center_shift_inverse_X==================================
     function: shift back positions of N (num) particles (p) based on chain center-of-mass
     p: particle list (read data and return shifted list);
  */
  void center_shift_inverse_X(particle *p) {
    for (int i=0;i<num;i++) {
      p[i].pos[0] += cm.pos[0];
      p[i].pos[1] += cm.pos[1];
      p[i].pos[2] += cm.pos[2];
    }
  }

  /* center_shift_X==================================
     function: shift positions and velocities of N (num) particles (p) based on chain center-of-mass
     p: particle list (read data and return shifted list);
  */
  void center_shift_X(particle *p) {
    for (int i=0;i<num;i++) {
      p[i].pos[0] -= cm.pos[0];
      p[i].pos[1] -= cm.pos[1];
      p[i].pos[2] -= cm.pos[2];
    }
  }

  /* set_abg
     function: set time step integration parameter alpha, beta and gamma
     argument: a: alpha
               b: beta
               g: gamma
               tmin: minimum time step
  */
  void set_abg(const double a, const double b, const double g, const double tmin) {
    alpha = a;
    beta = b;
    gamma = g;
    dtmin = tmin;
  }


  /* Leapfrog_step_forward
     function: integration n steps
     argument: s: whole step size
               n: number of step division, real step size = (s/n), for Leapfrog integration, it is X(s/2n)V(s/n)X(s/n)V(s/n)..X(s/2n)
               p: particle list (the position and velocity will be updated with center-mass-corrected)
               np: number of perturber particle lists
               pext: perturber particle lists
               force: external force (not perturber forces which are calculated in pert_force)
               check_flag: true: check link every step; false: check link at then end
///               force: external force array
///               upforce: void (const particle * p, const particle *pext, double3* force). function to calculate force based on p and pext, return to force
  */             
  //  void Leapfrog_step_forward(const double s, const int n, particle* p, particle* pext, double3* force, ext_force<particle> upforce) {
  void Leapfrog_step_forward(const double s, const int n, particle* p, const int np=0, const particle* pext=NULL, const double3* force=NULL, bool check_flag=false) {
    double ds = s/double(n);
    double3* ave_v=new double3[num];  // average velocity
    bool fpf = false; // perturber force indicator
    if (np>0) fpf = true;

    for (int i=0;i<n;i++) {
      // half step forward for X, t (dependence: Ekin, B, w, V)
      step_forward_X(ds/2.0);

      // perturber force
      if (fpf) {
        // resolve X to p.v (dependence: X, cm.mass)
        resolve_X(p);
        // get original position first, update p.x (dependence: X, cm.x)
        center_shift_inverse_X(p);
        // Update perturber force pf (dependence: pext, original-frame p.x, p.mass)
        pert_force(p,np,pext);
        // reset position to center-of-mass frame, update p.x 
        center_shift_X(p);
      }
      
      // Update rjk, A, Pot, dmdr, W for half X (dependence: pf, force, p.mass, X, Wjk)
      calc_rAPW(p,force); 

      // update chain list order if necessary, update list, X, V (dependence: rjk, V)
      if (num>2&&check_flag) update_link();

      // Step forward V and get time step dt(V)
      double dvt = step_forward_V(ds);
      
      // Get averaged velocity, update p.x, p.v, ave_v (dependence: X, V)
      resolve_XV(p,ave_v);

      // forward B and w (dependence: dt(V), ave_v, p.mass, p.v, force, dmdr, pf)
      step_forward_Bw(dvt,ave_v,p,force,fpf);

      // Calcuale Kinetic energy (dependence: p.mass, p.v)
      calc_Ekin(p);

      // half step forward for X (dependence: Ekin, B, w, V)
      step_forward_X(ds/2.0);

    }

   // resolve X at last, update p.x (dependence: X)
    resolve_X(p);

    // Update rjk, A, Pot, dmdr, W (notice A will be incorrect since pf is not updated)
    calc_rAPW(p,force); 
    
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
    if(num>2&&!check_flag)  update_link();
    
//    else if (num==2) {
//      calc_Pot(p); // for energy check (dependence: p.mass, p.x)
//    }
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
    if (dt2==0)
      if (dt1==0) return tik1;
      else {
        std::cerr<<"Error!: T_i,k-1 - T_i-1,k-2 = 0 but T_i,k-1 - T_i-1,k-1 = "<<dt1<<" (T_i,k-1 = "<<tik1<<"; T_i-1,k-1 = "<<ti1k1<<"; h_i-k/h_i = "<<hr<<")"<<std::endl;
        abort();
      }
    else return tik1 + (tik1 - ti1k1)/(hr*hr * (1 - dt1/dt2) - 1);
  }
  
  /* extrapolation_integration
     function: extrapolation method to get accurate integration
     argument: s: whole step size
               p: particle list (the position and velocity will be updated with center-mass-corrected)
               error: relative error requirement for extrapolation
               itermax: maximum times for iteration.
               methods: 1: Romberg method; others: Rational interpolation method
               np: number of perturber particle lists
               pext: perturber particle lists
               force: external force (not perturber forces which are calculated in pert_force)
               check_flag: true: check link every step; false: check link at then end
     return:   iteration count (start from 1, which means no iteration)
   */
  int extrapolation_integration(const double s, particle* p, const double error=1E-8, const int itermax=10, const int methods=1, const int np=0, const particle* pext=NULL, const double3* force=NULL, bool check_flag=false) {
    // for storage
    double   t0,ttemp,t1[itermax];
    double   B0,Btemp,B1[itermax];
    double   w0,wtemp,w1[itermax];
    double   Ekin0;
    //particle* p0 = new particle[num];

    double3* X1[itermax];
    double3* V1[itermax];
    for (int i=0;i<itermax;i++) {
      X1[i] = new double3[num];
      V1[i] = new double3[num];
    }
    double3* X0 = new double3[num];
    double3* V0 = new double3[num];
    double3* Xtemp = new double3[num];
    double3* Vtemp = new double3[num];

    double3 CX;
    double cxerr=error+1.0;

    // backup initial values
    t0 = t;
    B0 = B;
    w0 = w;
    Ekin0 = Ekin;
    memcpy(X0,X,3*num*sizeof(double));
    memcpy(V0,V,3*num*sizeof(double));
    //memcpy(p0,p,num*sizeof(particle));

    // first step
    Leapfrog_step_forward(s,1,p,np,pext,force,false);

    // relative position vector between first and last particle for phase error check
    int k0 = list[0];
    int kn = list[num-1];
    memcpy(CX,rjk[k0][kn],3*sizeof(double));

    int intcount = 0; // iteration counter
    int step = 1; //substep size
    while (cxerr > error) {
      if (intcount == 0) {
        // storage the results
        t1[0] = t;
        B1[0] = B;
        w1[0] = w;
        memcpy(X1[0],X,3*num*sizeof(double));
        memcpy(V1[0],V,3*num*sizeof(double));
      }
      intcount++;
      if (intcount == itermax) {
        std::cerr<<"Error: maximum iteration step number "<<itermax<<" reached, but phase error "<<cxerr<<" is larger than criterion "<<error<<std::endl;
        abort();
      }
      step = 2*step;

      // reset the initial data
      t = t0;
      B = B0;
      w = w0;
      Ekin = Ekin0;
      memcpy(X,X0,3*num*sizeof(double));
      memcpy(V,V0,3*num*sizeof(double));

      Leapfrog_step_forward(s,step,p,np,pext,force,false);

      if (methods==1) {
        // Using Romberg method
        // first step extrapolate, storage the temp data in last index of array
        /*
          T_n-1,0 [0]
                     -> T_n,1 [n]
          T_n,0 [Chain]
        */
        // H_i-k/H_i
        int hr = 2;
        
        t1[intcount] = Romberg_recursion_formula(t1[0],t,hr);
        B1[intcount] = Romberg_recursion_formula(B1[0],B,hr);
        w1[intcount] = Romberg_recursion_formula(w1[0],w,hr);
        for (int i=0;i<num;i++) {
          for (int k=0;k<3;k++) {
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
          memcpy(X1[0],X,3*num*sizeof(double));
          memcpy(V1[0],V,3*num*sizeof(double));

          // iteration to get final result
          for (int j=1; j<intcount; j++) {
            //templately storage new results to ttemp
            /*
              T_n-1,j [j]
                       -> T_n,j+1 [temp]
              T_n,j [n]
             */
            hr = 2*hr;
            ttemp = Romberg_recursion_formula(t1[j],t1[intcount],hr);
            Btemp = Romberg_recursion_formula(B1[j],B1[intcount],hr);
            wtemp = Romberg_recursion_formula(w1[j],w1[intcount],hr);
            for (int i=0;i<num;i++) {
              for (int k=0;k<3;k++) {
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
            memcpy(X1[j],X1[intcount],3*num*sizeof(double));
            memcpy(V1[j],V1[intcount],3*num*sizeof(double));
   
            //shift temp data to index = intcount.
            /*
              [n] << T_n,j+1 [temp]
             */
            t1[intcount] = ttemp;
            B1[intcount] = Btemp;
            w1[intcount] = wtemp;
            memcpy(X1[intcount],Xtemp,3*num*sizeof(double));
            memcpy(V1[intcount],Vtemp,3*num*sizeof(double));
          }
        }
      }
      else {
        // Using Rational interpolation method
        // additional template storage
        double tt1, Bt1, wt1;
        double3* Xt1 = new double3[num];
        double3* Vt1 = new double3[num];
        int hr = 2;

        /*
           T_n-1,0 [0]
        0               -> T_n,1 [n]
           T_n,0 [Chain]
        */
        t1[intcount] = Rational_recursion_formula(0,t1[0],t,hr);
        B1[intcount] = Rational_recursion_formula(0,B1[0],B,hr);
        w1[intcount] = Rational_recursion_formula(0,w1[0],w,hr);
        for (int i=0;i<num;i++) {
          for (int k=0;k<3;k++) {
            X1[intcount][i][k] = Rational_recursion_formula(0,X1[0][i][k],X[i][k],hr);
            V1[intcount][i][k] = Rational_recursion_formula(0,V1[0][i][k],V[i][k],hr);
          }
        }

        // if step number > 1, need iteration
        if (intcount>1) {
          // iteration to get final result
          for (int j=1; j<intcount; j++) {
            //templately storage new results to ttemp
            /*
                           T_n-1,j [j]
             T_n-1,j-1 [j-1]          -> T_n,j+1 [temp]
                           T_n,j [n]
             */
            hr = 2*hr;
            ttemp = Rational_recursion_formula(t1[j-1],t1[j],t1[intcount],hr);
            Btemp = Rational_recursion_formula(B1[j-1],B1[j],B1[intcount],hr);
            wtemp = Rational_recursion_formula(w1[j-1],w1[j],w1[intcount],hr);
            for (int i=0;i<num;i++) {
              for (int k=0;k<3;k++) {
                Xtemp[i][k] = Rational_recursion_formula(X1[j-1][i][k],X1[j][i][k],X1[intcount][i][k],hr);
                Vtemp[i][k] = Rational_recursion_formula(V1[j-1][i][k],V1[j][i][k],V1[intcount][i][k],hr);
              }
            }

            if (j==1) {
              // update zero order
              // [0] << T_n,0
              t1[0] = t;
              B1[0] = B;
              w1[0] = w;
              memcpy(X1[0],X,3*num*sizeof(double));
              memcpy(V1[0],V,3*num*sizeof(double));
            }
            else {
              // update j-1 order
              // [j-1] << T_n,j-1 [t1]
              t1[j-1] = tt1;
              B1[j-1] = Bt1;
              w1[j-1] = wt1;
              memcpy(X1[j-1],Xt1,3*num*sizeof(double));
              memcpy(V1[j-1],Vt1,3*num*sizeof(double));
            }

            //storage previous extrapolation data in template position t1
            // [t1] << T_n,j [n]
            tt1 = t1[intcount];
            Bt1 = B1[intcount];
            wt1 = w1[intcount];
            memcpy(Xt1,X1[intcount],3*num*sizeof(double));
            memcpy(Vt1,V1[intcount],3*num*sizeof(double));
   
            //shift temp data to index = intcount.
            // [n] << T_n,j+1 [temp]
            t1[intcount] = ttemp;
            B1[intcount] = Btemp;
            w1[intcount] = wtemp;
            memcpy(X1[intcount],Xtemp,3*num*sizeof(double));
            memcpy(V1[intcount],Vtemp,3*num*sizeof(double));
          }
        }
      }

      // set final results back to chain array
      t = t1[intcount];
      B = B1[intcount];
      w = w1[intcount];
      memcpy(X,X1[intcount],3*num*sizeof(double));
      memcpy(V,V1[intcount],3*num*sizeof(double));

      // resolve particle
      resolve_XV(p);
      // recalculate the energy
      calc_Ekin(p);
      //calc_Pot(p);
      calc_rAPW(p);

      // phase error calculation
      double* r0n = rjk[k0][kn];
      double dcx1 = r0n[0] - CX[0];
      double dcx2 = r0n[1] - CX[1];
      double dcx3 = r0n[2] - CX[2];
      cxerr = std::sqrt(dcx1*dcx1 + dcx2*dcx2 + dcx3*dcx3)/r0n[3];
      memcpy(CX,r0n,3*sizeof(double));
#ifdef DEBUG
      std::cerr<<std::setprecision(14)<<"Iteration: "<<intcount<<" error = "<<cxerr;
      std::cerr<<" energy error "<<Ekin-Pot+B<<std::endl;
#endif 
    }
    
    return intcount+1;
  }

  /* GetTime
     function: returen physical time
     return: t
  */
  double getTime() {
    return t;
  }
  /* get_Ekin
     function: get kinetic energy
     return: Ekin
  */
  double getEkin() {
    return Ekin;
  }

  /* get_Pot
     function: get potetnial energy (negative value)
     return: Pot
  */
  double getPot() {
    return -Pot;
  }

  /* get_B
     function: get potetnial energy (negative value)
     return: Pot
  */
  double getB() {
    return B;
  }
  
  /* get_w
     function: get potetnial energy (negative value)
     return: Pot
  */
  double getw() {
    return w;
  }

  /* get_W
     function: get potetnial energy (negative value)
     return: Pot
  */
  double getW() {
    return W;
  }
  

#ifdef DEBUG
  /* set_X
     function: modify X
     argument: i: the index of X
               k: axis-0:x,1:y,2:z
               value: new value
   */
  void set_X(const int i, const int k, const double value) {
    X[i][k] = value;
  }


  //  template <class particle>
  void update_rAPW(const particle* p, const double3* force) {
    calc_rAPW(p,force);
  }
  
  /* print
     function: print all data
     argument: width of each output value
  */
  //  template <class particle>
  void print(const int width) {
    if (width<=0) {
      std::cerr<<"Error: width should be larger than zero!\n";
      abort();
    }
    char xyz[4]={'x','y','z','r'};
    std::cout<<"---- particle list------\n ";
    for (int i=0;i<num;i++) std::cout<<std::setw(width)<<list[i];
    std::cout<<"\n----- relative position X ------\n";
    for (int k=0;k<3;k++) {
      std::cout<<xyz[k];
      for (int i=0;i<num-1;i++) std::cout<<std::setw(width)<<X[i][k];
      std::cout<<std::endl;
    }
    std::cout<<"\n----- relative velocity V ------\n";
    for (int k=0;k<3;k++) {
      std::cout<<xyz[k];
      for (int i=0;i<num-1;i++) std::cout<<std::setw(width)<<V[i][k];
      std::cout<<std::endl;
    }
    std::cout<<"\n----- mass coefficients Matrix Wjk -----\n";
    for (int i=0;i<num;i++){
      std::cout<<' ';
      for (int j=0;j<num;j++) std::cout<<std::setw(width)<<Wjk[i][j];
      std::cout<<std::endl;
    }
    std::cout<<"\n----- relative distances Matrix rjk -----\n";
    for (int k=0;k<4;k++){
      std::cout<<"----------"<<xyz[k]<<"----------\n";
      for (int i=0;i<num;i++){
        std::cout<<' ';
        for (int j=0;j<num;j++) std::cout<<std::setw(width)<<rjk[i][j][k];
        std::cout<<std::endl;
      }
    }
    std::cout<<"\n----- Acceleration A ------\n";
    for (int k=0;k<3;k++) {
      std::cout<<xyz[k];
      for (int i=0;i<num;i++) std::cout<<std::setw(width)<<acc[i][k];
      std::cout<<std::endl;
    }
    std::cout<<"\n----- part omega / part rk ------\n";
    for (int k=0;k<3;k++) {
      std::cout<<xyz[k];
      for (int i=0;i<num;i++) std::cout<<std::setw(width)<<dmdr[i][k];
      std::cout<<std::endl;
    }
    std::cout<<"\n----- system parameters ------\n"
             <<"---Center-of-mass data: \n"
             <<"mass:"<<std::setw(width)<<cm.mass<<std::endl
             <<"(particle*)cm.pos:"
             <<std::setw(width)<<cm.pos[0]
             <<std::setw(width)<<cm.pos[1]
             <<std::setw(width)<<cm.pos[2]<<std::endl
             <<"(particle*)cm.vel:"
             <<std::setw(width)<<cm.vel[0]
             <<std::setw(width)<<cm.vel[1]
             <<std::setw(width)<<cm.vel[2]<<std::endl
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

