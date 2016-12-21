#pragma once

#include <cmath>

namespace TB {

  //! Calculate parameters of two-body motion
  /*! Calcuate the basic parameters of two-body motions
    @param[in] semi: semi-major axis
    @param[in] peri: period
    @param[in] ecc:  eccentricity
    @param[in] m: total mass of two particles
    @param[in] dx: relative position [3]
    @param[in] dv: relative velocity [3]
   */
  void calc_pars_two(double &semi, double &peri, double &ecc, const double m, const double dx[3], const double dv[3]) {
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

}
