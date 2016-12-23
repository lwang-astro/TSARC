#pragma once

#include <cstdlib>
#include <cmath>

namespace NTA {

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

  //! Calculate parameters of two-body motion
  /*! Calcuate the basic parameters of two-body motions
    @param[in] semi: semi-major axis
    @param[in] peri: period
    @param[in] ecc:  eccentricity
    @param[in] m: total mass of two particles
    @param[in] dx: relative position [3]
    @param[in] dv: relative velocity [3]
   */
  void calc_kepler_orbit_par(double &semi, double &peri, double &ecc, const double m, const double dx[3], const double dv[3]) {
    const double dr2  = dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2];
    const double dv2  = dv[0]*dv[0] + dv[1]*dv[1] + dv[2]*dv[2];
    const double rdot = dx[0]*dv[0] + dx[1]*dv[1] + dx[2]*dv[2];
    const double dr   = std::sqrt(dr2);

    semi = 1.0/(2.0/dr - dv2/m);

    const double dr_semi = 1.0 - dr/semi;
    ecc  = std::sqrt(dr_semi*dr_semi + rdot*rdot/(m*semi));

    const double twopi = 4.0*std::atan(1.0);
    peri = twopi*std::abs(semi)*std::sqrt(std::abs(semi)/m);
#ifdef DEBUG    
    std::cerr<<"m="<<m<<" dx="<<dx[0]<<" "<<dx[1]<<" "<<dx[2]<<" dv="<<dv[0]<<" "<<dv[1]<<" "<<dv[2]<<" semi="<<semi<<" ecc="<<ecc<<" peri="<<peri<<std::endl;
#endif    
  }

  //! Newtonian acceleration and \f$\partial W_{ij}/\partial \mathbf{x}_i\f$ from particle j to particle i (function type of \link #ARC::pair_AW \endlink) 
  /*!          @param[out] Aij: Newtonian acceleration vector for i particle from particle j. \f$Aij[1:3] = m_i m_j xij[1:3] / |xij|^3 \f$.
    @param[out] Pij: Newtonian potential of i from j. \f$ Pij = -m_i m_j /|xij| \f$
    @param[out] pWij: TTL time transformation function partial derivates (component from j to i) \f$\partial W_{ij}/\partial \mathbf{x}_i\f$ (used for TTL method). \f$pWij[1:3] = mm_{ij} xij[1:3] /|xij|^3 \f$. (Total value is \f$\frac{\partial W}{\partial \mathbf{x}_i} = \sum_{j} mm_{ij} \mathbf{x}_{ij}/|\mathbf{x}_{ij}|^3\f$)
    @param[out] Wij: TTL time transformation function component with i,j (used for TTL method) \f$Wij = mm_{ij} /|xij|^3\f$ total value is \f$ W = \sum_{i<j} mm_{ij} /|xij| \f$
    @param[in] xij: relative position vector [1:3] \f$ \mathbf{x_j} - \mathbf{x_i} \f$
    @param[in] mi: particle i mass.
    @param[in] mj: particle j mass.
    @param[in] pars: Newtionian_pars type data with members:
    - mm2: smooth mass coefficient \f$ \sum_{i<j} m_i m_j /(N(N-1)/2) \f$ (can be calculated by calc_mm2()); \n
    - epi: adiustable parameter. 
    1) If epi>0: if \f$m_i m_j < epi mm2\f$:  \f$ mm_{ij} = mm2\f$  else: \f$ mm_ij = 0\f$
    2) If epi<0: \f$mm_{ij} = m_i m_j\f$.\n
  */
  void Newtonian_AW (double Aij[3], double &Pij, double pWij[3], double &Wij, const double xij[3], const double &mi, const double &mj, const Newtonian_pars* pars) {

    // safetey check
    if (pars==NULL) {
      std::cerr<<"Error: Newtonian_pars is not set in chain!\n";
      abort();
    }
    
    // distance
    double rij = std::sqrt(xij[0]*xij[0]+xij[1]*xij[1]+xij[2]*xij[2]);  

    // smooth coefficients
    double mm2=pars->mm2;
    double epi=pars->epi;

    // mass parameters
    double mimj = mi*mj; // m_i*m_i
    double mmij;
    if (mm2>0 && epi>0) {
      // Wij = mm2 if m_i*m_i < epi*m'^2; 0 otherwise;
      if (mimj<epi*mm2) mmij = mm2;
      else mmij = 0;
    }
    else {
      mmij = mimj;    // Wij = m_i*m_i
    }
  
    Pij = - mimj / rij;  // Potential energy
    Wij = mmij / rij;   // Transformation coefficient
        
    // Acceleration
    double rij3 = rij*rij*rij;
    double mor3 = mj / rij3;
    Aij[0] = mor3 * xij[0];
    Aij[1] = mor3 * xij[1];
    Aij[2] = mor3 * xij[2];

    // dW/dr
    mor3 = mmij / rij3;
    pWij[0] = mor3 * xij[0];
    pWij[1] = mor3 * xij[1];
    pWij[2] = mor3 * xij[2];
  
  }

  //! Newtonian acceleration from particle p to particle i (function type of ::ARC::pair_Ap)
  /*! 
    @param[out]  Aij: acceleration vector. \f$Aij[1:3] = m_i m_p (xp[1:3]-xi[1:3]) / |xp-xi|^3 \f$.
    @param[out]  Pij: potential. \f$ Pij = - m_i m_p /|xp-xi|^3\f$
    @param[in]  xi: position vector i.
    @param[in]  xp: position vector p.
    @param[in]  mi: particle mass i.
    @param[in]  mp: particle mass p.
    @param[in]  pars: Newtonian_pars type data (not used)
  */
  void Newtonian_Ap (double Aij[3], double &Pij, const double xi[3], const double xp[3], const double &mi, const double &mp, const Newtonian_pars* pars){
    double dx = xp[0] - xi[0];
    double dy = xp[1] - xi[1];
    double dz = xp[2] - xi[2];

    double dr2 = dx*dx + dy*dy + dz*dz;
    double dr  = std::sqrt(dr2);
    double dr3 = dr*dr2;

    Aij[0] = mp * dx / dr3;
    Aij[1] = mp * dy / dr3;
    Aij[2] = mp * dz / dr3;

    Pij = - mi*mp / dr;
  
  }

  //! Newtonian two-body kepler period
  /*! Use calc_kepler_orbit_par() to obtain the two-body kepler period
    @param[in] m1: mass of particle 1
    @param[in] m2: mass of particle 2
    @param[in] X: relative position vector
    @param[in] V: relative velocity vector
    @param[in] pars: Newtonian_pars type data (not used)
    \return Return the period
  */
  double Newtonian_kepler_period (const double m1, const double m2, const double X[3], const double V[3], const Newtonian_pars* pars) {
    double semi, ecc, peri;
    calc_kepler_orbit_par(semi, peri, ecc, m1+m2, X, V);
    return peri;
  }

}
