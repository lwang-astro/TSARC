
#include "AR.h"
#include <initial.h>
#include <iostream>
#include <fstream>

int main(int argc, char **argv){
  pars_initial init(".debug_confg");
  init.add("n","particle number",(int)3);
  init.add("filename","particle data filename (m,x,v)","input");
  init.add("width","print width",(int)15);
  init.add("nstep","total step number",(int)3000);
  init.add("stepsize","step size",(double)(std::pow(0.5,6)));
  init.add("error","error limit",(double)1E-7);
  init.add("itermax","Maximum iteration steps",(int)10);
  init.add("stepdiv","sub-step number",(int)128);
  init.add("methods","1: Extrapolation; 2: LF",(int)2);
  

  init.initial(argc,argv);
  int n=init.geti("n");
  std::fstream fs;
  fs.open(init.gets("filename").c_str(),std::fstream::in);
  int w=init.geti("width");

  chain<Particle> c(n);
  Particle *p=new Particle[n];
  for (int i=0;i<n;i++) {
    double x,y,z,vx,vy,vz,m;
    fs>>m>>x>>y>>z>>vx>>vy>>vz;
    p[i]=Particle(m,x,y,z,vx,vy,vz);
  }

  char xyz[3]={'x','y','z'};

  if (false) {
    std::cout<<"-------Particle:\n"
             <<"---------mass:-------\n "
             <<std::setw(w)<<p[0].mass
             <<std::setw(w)<<p[1].mass
             <<std::setw(w)<<p[2].mass;
    std::cout<<"\n---------Position:-----------\n";
    for (int i=0;i<3;i++) {
      std::cout<<xyz[i];
      for (int j=0;j<3;j++) {
        std::cout<<std::setw(w)<<p[j].pos[i];
      }
      std::cout<<std::endl;
    }
    std::cout<<"---------Velocity:-----------\n";
    for (int i=0;i<3;i++) {
      std::cout<<xyz[i];
      for (int j=0;j<3;j++) {
        std::cout<<std::setw(w)<<p[j].pos[i];
      }
      std::cout<<std::endl;
    }
  }

  c.set_abg(1.0,0.0,0.0,std::pow(0.5,40));
  c.init(0.0,n,p);
  //  std::cout<<std::setprecision(w-2);
  //  c.print(w);
  int nstep=init.geti("nstep");
  int nsubstep=init.geti("stepdiv");
  int itermax=init.geti("itermax");
  int sw=init.geti("methods");
  double err=init.getd("error");
  double s=init.getd("stepsize");
  int stepcount[32]={};
  for (int i=0;i<nstep;i++) {
    std::cout<<c.getTime();
    for (int j=0;j<n;j++) {
      for (int k=0;k<3;k++) {
        std::cout<<std::setw(w)<<p[j].pos[k];
      }
    }
    std::cout<<std::setw(w)<<c.getEkin()+c.getPot()+c.getB()<<std::setw(w)<<c.getEkin()<<std::setw(w)<<c.getPot()
    <<std::setw(w)<<c.getB()<<std::setw(w)<<c.getw()<<std::endl;
    if (sw==2) c.Leapfrog_step_forward(s,nsubstep,p);
    else stepcount[c.extrapolation_integration(s,p,err,itermax)-1]++;
  }

  if (sw!=2) {
    std::cerr<<"Step histogram:\n I\tCount\n";
    for (int i=0;i<32;i++) {
      std::cerr<<i+1<<"\t"<<stepcount[i]<<std::endl;
    }
  }
  //c.print(w);

//  c.set_X(0,1,2);
//  c.set_X(1,1,-2);
//  c.update_rAPW(p,f);
//  c.update_link();
//  c.print(w);

  return 0;
}
  
