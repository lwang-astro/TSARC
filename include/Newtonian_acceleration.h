#pragma once

#include <cstdlib>
#include <cmath>
#include "particle.h"

//! Namespace for Newtonian Interaction related functions
/*!
  In the special applications of N-body simulations with Newtonian gravity, the Newtonian acceleration functions should be defined. 
  The corresponding timescale to determing the next integrating step should also be provided.
  This is the namespace where all these functions are provided.
*/
namespace NTA {

  typedef double double3[3];

  //! Newtonian Interaction and time transformation function parameter class
  /*!
    This class store the smooth mass coefficients for TTL methods
   */
  class Newtonian_pars{
  public:
    double mm2; /// smooth particle coefficient
    double epi; /// adjustment parameter

    //! constructor
    /*! #mm2 is set to zero and #epi is set to -1
     */
    Newtonian_pars(): mm2(0.0), epi(-1.0) {}

    //! Calculate smooth particle coefficient
    /*! Get smooth particle coefficient \f$ mm2 = \sum m_i  m_j /(N (N-1)/2) (i<j) \f$ where \f$m_i\f$ and \f$m_j\f$ are particle masses
      @param[in] mass: particle mass array
      @param[in] n: number of particle
      \return smooth particle coefficient
    */
    void calc_mm2(const double mass[], const std::size_t n) {
      // calcualte m'^2
      mm2 = 0;
      for (std::size_t i=0;i<n;i++) {
        for (std::size_t j=i+1;j<n;j++) {
          mm2 += mass[i]* mass[j];
        }
      }
      mm2 /= n * (n - 1.0) / 2.0;
    }    

  };

  //! Newtonian acceleration and \f$\partial W_{ij}/\partial \mathbf{x}_i\f$ from particle j to particle i (function type of \link #ARC::pair_AW \endlink) 
  /*!          @param[out] Aij: Newtonian acceleration vector for i particle from particle j. \f$Aij[1:3] = m_i m_j xij[1:3] / |xij|^3 \f$.
    @param[out] Pij: Newtonian potential of i from j. \f$ Pij = -m_i m_j /|xij| \f$
    @param[out] pWij: TTL time transformation function partial derivates (component from j to i) \f$\partial W_{ij}/\partial \mathbf{x}_i\f$ (used for TTL method). \f$pWij[1:3] = mm_{ij} xij[1:3] /|xij|^3 \f$. (Total value is \f$\frac{\partial W}{\partial \mathbf{x}_i} = \sum_{j} mm_{ij} \mathbf{x}_{ij}/|\mathbf{x}_{ij}|^3\f$)
    @param[out] Wij: TTL time transformation function component with i,j (used for TTL method) \f$Wij = mm_{ij} /|xij|^3\f$ total value is \f$ W = \sum_{i<j} mm_{ij} /|xij| \f$
    @param[in] xij: relative position vector [1:3] \f$ \mathbf{x_j} - \mathbf{x_i} \f$
    @param[in] pi: particle i (get mass)
    @param[in] pj: particle j (get mass)
    @param[in] pars: Newtionian_pars type data with members:
    - mm2: smooth mass coefficient \f$ \sum_{i<j} m_i m_j /(N(N-1)/2) \f$ (can be calculated by calc_mm2()); \n
    - epi: adiustable parameter. 
    1) If epi>0: if \f$m_i m_j < epi mm2\f$:  \f$ mm_{ij} = mm2\f$  else: \f$ mm_ij = 0\f$
    2) If epi<0: \f$mm_{ij} = m_i m_j\f$.\n
    \return status: 0 for normal cases; 1 for the case when two particles have same positions
  */
  int Newtonian_AW (double Aij[3], double &Pij, double pWij[3], double &Wij, const double xij[3], const Particle &pi, const Particle &pj, Newtonian_pars* pars) {

    // safetey check
    if (pars==NULL) {
      std::cerr<<"Error: Newtonian_pars is not set in chain!\n";
      abort();
    }
    
    // distance
    double rij = std::sqrt(xij[0]*xij[0]+xij[1]*xij[1]+xij[2]*xij[2]);

    if (rij==0) return 1;

    // smooth coefficients
    double mm2=pars->mm2;
    double epi=pars->epi;

    // mass parameters
    double mimj = pi.getMass()*pj.getMass(); // m_i*m_j
    double mmij;
    if (mm2>0 && epi>0) {
      // Wij = mm2 if m_i*m_j < epi*m'^2; 0 otherwise;
      if (mimj<epi*mm2) mmij = mm2;
      else mmij = 0;
    }
    else {
      mmij = mimj;    // Wij = m_i*m_j
    }
  
    Pij = - mimj / rij;  // Potential energy
    Wij = mmij / rij;   // Transformation coefficient
        
    // Acceleration
    double rij3 = rij*rij*rij;
    double mor3 = pj.getMass() / rij3;
    Aij[0] = mor3 * xij[0];
    Aij[1] = mor3 * xij[1];
    Aij[2] = mor3 * xij[2];

    // dW/dr
    mor3 = mmij / rij3;
    pWij[0] = mor3 * xij[0];
    pWij[1] = mor3 * xij[1];
    pWij[2] = mor3 * xij[2];

    return 0;
  }

  //! Newtonian acceleration from perturber pert to particle p (function type of ::ARC::ext_Acc)
  /*! 
    @param[out] acc: acceleration vector from pert to p. \f$ \sum_j Aij[1:3] = m_i m_p (xp[1:3]-xi[1:3]) / |xp-xi|^3 \f$.
    @param[in] t: time step for prediction of pert particles
    @param[in] p: particle array
    @param[in] np: number of particles
    @param[in] pert: perturber particle array
    @param[in] pertf: perturrber force array for prediction
    @param[in] npert: number of perturbers
    @param[in] pars: extra parameters (not used)
  */
  template<class particle>  
  void Newtonian_extAcc(double3 *acc, const double t, particle *p, const int np, particle *pert, double3 *pertf, const int npert, Newtonian_pars *pars) {
      double3 xp[npert];
      for (int i=0; i<npert; i++) {
          const double* r = pert[i].getPos();
          const double* v = pert[i].getVel();
          double dt2 = t*t;
          for(int j=0; j<3; j++) {
              xp[i][j] = r[j]+v[j]*t + 0.5*pertf[i][j]*dt2;
          }
      }
      
      for (int i=0; i<np; i++) {
          const double* xi = p[i].getPos();
          acc[i][0] = acc[i][1] = acc[i][2] = 0.0;
          for (int j=0; j<npert; j++) {
      
              double mp = pert[i].getMass();
              double dx = xp[i][0] - xi[0];
              double dy = xp[i][1] - xi[1];
              double dz = xp[i][2] - xi[2];

              double dr2 = dx*dx + dy*dy + dz*dz;
              double dr  = std::sqrt(dr2);
              double dr3 = dr*dr2;

              acc[i][0] = mp * dx / dr3;
              acc[i][1] = mp * dy / dr3;
              acc[i][2] = mp * dz / dr3;

          }
      }
              //Pij = - pi.getMass()*mp / dr;
  }

  //! Newtonian two-body kepler period
  /*! Use calc_kepler_orbit_par() to obtain the two-body kepler period
    @param[in] m1: mass of particle 1
    @param[in] m2: mass of particle 2
    @param[in] dx: relative position vector
    @param[in] dv: relative velocity vector
    @param[in] pars: Newtonian_pars type data (not used)
    \return If the orbit is close, return the period, otherwise return the approximately 10% of free-fall time.
  */
  double Newtonian_kepler_period (const double m1, const double m2, const double dx[3], const double dv[3], Newtonian_pars* pars) {
    double semi;
    const double dr2  = dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2];
    const double dv2  = dv[0]*dv[0] + dv[1]*dv[1] + dv[2]*dv[2];
    const double dr   = std::sqrt(dr2);
    const double m    = m1+m2;

    semi = 1.0/(2.0/dr - dv2/m);

    if (semi<0) {
      const double peri = 0.1*std::sqrt(dr2*dr/(2.0*m));
      //      std::cout<<"dr="<<dr<<" semi="<<semi<<" peri="<<peri<<std::endl;
      return peri;
    }
    else {
      const double rvdot = dx[0]*dv[0] + dx[1]*dv[1] + dx[2]*dv[2];
      const double dr_semi = 1.0 - dr/semi;
      const double ecc  = std::sqrt(dr_semi*dr_semi + rvdot*rvdot/(m*semi));

      const double twopi= 6.28;
      const double peri = twopi*std::abs(semi)*std::sqrt(std::abs(semi)/m);

      //      std::cout<<"dr="<<dr<<" semi="<<semi<<" ecc="<<ecc<<" peri="<<peri<<std::endl;
      return std::max(std::sqrt(std::abs(1.0-ecc)),0.01)*peri;
    }
  }


  //! Calculate parameters of two-body motion
  /*! Calcuate the basic parameters of two-body motions
    @param[out] semi: semi-major axis
    @param[out] peri: period
    @param[out] ecc:  eccentricity
    @param[out] angle: three rotational angles: inclination, z-axis rotation, orbital-plane rotation
    @param[out] true_anomaly: true anomaly
    @param[out] ecc_anomaly: eccentricity anomaly
    @param[out] mean_anomaly: mean anomaly
    @param[in] m: total mass of two particles
    @param[in] dx: relative position [3]
    @param[in] dv: relative velocity [3]
   */
  void calc_kepler_orbit_par(double &semi, double &peri, double &ecc, double angle[3], double &true_anomaly, double &ecc_anomaly, double &mean_anomaly, const double m, const double dx[3], const double dv[3]) {
    const double dr2  = dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2];
    const double dv2  = dv[0]*dv[0] + dv[1]*dv[1] + dv[2]*dv[2];
    const double rvdot = dx[0]*dv[0] + dx[1]*dv[1] + dx[2]*dv[2];
    const double dr   = std::sqrt(dr2);

    semi = 1.0/(2.0/dr - dv2/m);

    const double dr_semi = 1.0 - dr/semi;
    ecc  = std::sqrt(dr_semi*dr_semi + rvdot*rvdot/(m*semi));

    const double twopi = 8.0*std::atan(1.0);
    peri = twopi*std::abs(semi)*std::sqrt(std::abs(semi)/m);

    const double rvcross[3]={dx[1]*dv[2]-dx[2]*dv[1],
                             dx[2]*dv[0]-dx[0]*dv[2],
                             dx[0]*dv[1]-dx[1]*dv[0]};

    // inclination to x-y plane
    const double rvcxy = std::sqrt(rvcross[0]*rvcross[0]+rvcross[1]*rvcross[1]);
    angle[0] = std::atan2(rvcxy, rvcross[2]);
    // z-axis rotation (counter-clock)
    angle[1] = std::atan2(rvcross[0], -rvcross[1]);
    if (angle[0]==0.0) angle[1]=0.0;

    const double cosinc = std::cos(angle[0]);
    const double sininc = std::sin(angle[0]);

    const double costhe = std::cos(angle[1]);
    const double sinthe = std::sin(angle[1]);
    
    // orbital-plane projection
    const double dr_op[3]={dx[0]*costhe + dx[1]*sinthe,
                           (-dx[0]*sinthe + dx[1]*costhe)*cosinc + dx[2]*sininc,
                           0};
    const double dv_op[3]={dv[0]*costhe + dv[1]*sinthe,
                           (-dv[0]*sinthe + dv[1]*costhe)*cosinc + dv[2]*sininc,
                           0};
//    const double drv_op[3]={rvcross[0]*costhe + rvcross[1]*sinthe,
//                          -(rvcross[0]*sinthe + rvcross[1]*costhe)*cosinc + rvcross[2]*sininc,
//                           0};
    const double h = std::sqrt(rvcxy*rvcxy + rvcross[2]*rvcross[2]);
    const double ecc_cos =  h/m*dv_op[1] - dr_op[0]/dr;
    const double ecc_sin = -h/m*dv_op[0] - dr_op[1]/dr;
    //    std::cout<<std::setprecision(15)<<"\ncos "<<ecc_cos<<" sin "<<ecc_sin<<" h "<<h<<" m "<<m<<"\n dv_op "<<dv_op[0]<<" "<<dv_op[1]<<" "<<dv_op[0]*dv_op[0]+dv_op[1]*dv_op[1]<<" dr_op "<<dr_op[0]<<" "<<dr_op[1]<<" "<<dr_op[0]*dr_op[0]+dr_op[1]*dr_op[1]<<"\n dr2 "<<dr2<<" dv2 "<<dv2<<" costhe "<<costhe<<" sinthe "<<sinthe<<" dx[0] "<<dx[0]<<" dx[1] "<<dx[1]<<" "<<dx[0]*dx[0]+dx[1]*dx[1]<<" "<<costhe*costhe+sinthe*sinthe<<std::endl;
    //    abort();

    // orbital-plane rotation
    angle[2] = std::atan2(ecc_sin, ecc_cos);
    if (ecc==0.0) angle[2]=0.0;

    // true anomaly
    true_anomaly = std::atan2(dr_op[1], dr_op[0]) - angle[2];

    // ecc anomaly
    const double cos_ecca = dr*std::cos(true_anomaly)/ semi + ecc;
    const double sin_ecca = dr*std::sin(true_anomaly)/(semi*sqrt(1.0 - ecc*ecc));
    ecc_anomaly = std::atan2(sin_ecca,cos_ecca);

    // mean anomaly
    mean_anomaly = ecc_anomaly - ecc*sin_ecca;
    
  }

  void kepler_orbit_generator(double x1[3], double x2[3], double v1[3], double v2[3], const double m1, const double m2, const double semi, const double ecc, const double angle[3], const double ecc_anomaly) {
    double m = m1+m2;
    double sine = std::sin(ecc_anomaly);
    double cose = std::cos(ecc_anomaly);

    double eoff=std::sqrt(1.0-ecc*ecc);
    double coff=std::sqrt(m/(semi*semi*semi));

    double dr[3]={ semi*(cose-ecc),               semi*eoff*sine,                     0.0};
    double dv[3]={-semi*coff*sine/(1.0-ecc*cose), semi*coff*eoff*cose/(1.0-ecc*cose), 0.0};

    //rotation
    const double sininc = std::sin(angle[0]);
    const double cosinc = std::cos(angle[0]);
    const double sinthe = std::sin(angle[1]);
    const double costhe = std::cos(angle[1]);
    const double sinorb = std::sin(angle[2]);
    const double cosorb = std::cos(angle[2]);

    double rm[3][3]={
      { cosorb*costhe - sinorb*sinthe*cosinc,
       -sinorb*costhe - cosorb*sinthe*cosinc,
        sinthe*sininc },
      { cosorb*sinthe + sinorb*costhe*cosinc,
       -sinorb*sinthe + cosorb*costhe*cosinc,
       -costhe*sininc },
      { sinorb*sininc,
        cosorb*sininc,
        cosinc}
    };

    double ndr[3]={rm[0][0]*dr[0] + rm[0][1]*dr[1] + rm[0][2]*dr[2],
                   rm[1][0]*dr[0] + rm[1][1]*dr[1] + rm[1][2]*dr[2],
                   rm[2][0]*dr[0] + rm[2][1]*dr[1] + rm[2][2]*dr[2]};
    double ndv[3]={rm[0][0]*dv[0] + rm[0][1]*dv[1] + rm[0][2]*dv[2],
                   rm[1][0]*dv[0] + rm[1][1]*dv[1] + rm[1][2]*dv[2],
                   rm[2][0]*dv[0] + rm[2][1]*dv[1] + rm[2][2]*dv[2]};

    for (int i=0; i<3; i++) {
      x1[i] = -m2/m * ndr[i];
      x2[i] =  m1/m * ndr[i];
      v1[i] = -m2/m * ndv[i];
      v2[i] =  m1/m * ndv[i];
    }
  }
  
}
