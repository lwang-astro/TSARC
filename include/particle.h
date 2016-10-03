#pragma once
#include <cassert>
#include <string.h>

#ifndef NAN_CHECK
#define NAN_CHECK(val) assert((val) == (val));
#endif

//basic particle structure==========================//
class Particle{
  double pos[3], vel[3];
  double mass;

public:
  //instruction=======================================//
  Particle() {} 
  Particle(const double &m, const double r[3], const double v[3]) {
    set(m,r,v);
  }
  Particle(const double &m, const double &rx, const double &ry, const double &rz, const double &vx, const double &vy, const double &vz) {
    set(m,rx,ry,rz,vx,vy,vz);
  }
  Particle(const Particle &a) {
    set(a);
  }

  // Get mass (required)
  const double getMass() const{
    return mass;
  }
  
  // Get position (required)
  const double* getPos() const{
    return pos;
  }

  // Get velocity (required)
  const double* getVel() const{
    return vel;
  }
  
  //set data=========================================//
  void set(const double &m, const double &rx, const double &ry, const double &rz, const double &vx, const double &vy, const double &vz){
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

  void set(const double &m, const double r[3], const double v[3]) {
    NAN_CHECK(m);
    NAN_CHECK(r[0]);
    NAN_CHECK(r[1]);
    NAN_CHECK(r[2]);
    NAN_CHECK(v[0]);
    NAN_CHECK(v[1]);
    NAN_CHECK(v[2]);

    mass=m;
    memcpy(pos,r,3*sizeof(double));
    memcpy(vel,v,3*sizeof(double));
  }
    

  void set(const Particle &a){
    mass = a.getMass();
    memcpy(pos,a.getPos(),3*sizeof(double));
    memcpy(vel,a.getVel(),3*sizeof(double));
  }

  //set position (required)
  void setPos(const double &x, const double &y, const double &z) {
    NAN_CHECK(x);
    NAN_CHECK(y);
    NAN_CHECK(z);
    
    pos[0] = x;
    pos[1] = y;
    pos[2] = z;
  }

  //set position (required)
  void setPos(const double r[3]) {
    NAN_CHECK(r[0]);
    NAN_CHECK(r[1]);
    NAN_CHECK(r[2]);

    memcpy(pos,r,3*sizeof(double));
  }

  //set velocity (required)
  void setVel(const double &vx, const double &vy, const double &vz) {
    NAN_CHECK(vx);
    NAN_CHECK(vy);
    NAN_CHECK(vz);
    
    vel[0] = vx;
    vel[1] = vy;
    vel[2] = vz;
  }

  //set position (required)
  void setVel(const double v[3]) {
    NAN_CHECK(v[0]);
    NAN_CHECK(v[1]);
    NAN_CHECK(v[2]);

    memcpy(vel,v,3*sizeof(double));
  }

  //set mass (required)
  void setMass(const double &m) {
    NAN_CHECK(m);

    mass = m;
  }

  void clear(){
    memset(pos,0,3*sizeof(double));
    memset(vel,0,3*sizeof(double));
    mass = 0.0;
  }
  
};
