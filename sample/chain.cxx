#include "AR.h"
#include "particle.h"
#include <iostream>
#include <fstream>
#include <unistd.h>
#include <string.h>
#include <stdlib.h>
#include <iomanip>

int main(int argc, char **argv){
  int n=3; //particle number
  int w=18; //print width
  int pre=10; //print digital precision
  int nstep=1000; // total step size
  int nsubstep=128;  // sub-step number if direct LF method is used
  int itermax=15;  //iteration maximum number for extrapolation methods
  char* sw=NULL;        // use extrapolation method, 'linear' for Romberg method; 'rational'  for rational interpolation method
  char* sq=NULL;        // extrapolation sequence, 'even' for {h,h/2,h/4,h/8...}; 'bs' for {h,h/2,h/3,h/4,h/6,h/8...}
  char* method=NULL;   // regularization methods, 'logh': Logarithmic Hamitonian; 'ttl': Time-transformed Leapfrog\n (logh)
  double err=1e-8; // phase error requirement
  double s=0.5;    // step size
  double dtmin=5.4e-20; // mimimum physical time step
  int copt;

  while ((copt = getopt(argc, argv, "a:d:N:w:n:s:k:i:m:M:e:p:h")) != -1)
    switch (copt) {
    case 'N':
      n = atoi(optarg);
      break;
    case 'w':
      w = atof(optarg);
      break;
    case 'n':
      nstep = atoi(optarg);
      break;
    case 's':
      s = atof(optarg);
      break;
    case 'k':
      nsubstep = atoi(optarg);
      break;
    case 'i':
      itermax = atoi(optarg);
      break;
    case 'm':
      sw = optarg;
      if (strcmp(sw,"linear")&&strcmp(sw,"rational")) {
        std::cerr<<"Extrapolation method "<<sw<<" not found!\n";
        abort();
      }
      break;
    case 'M':
      sq = optarg;
      if (strcmp(sq,"even")&&strcmp(sq,"bs")) {
        std::cerr<<"Extrapolation sequences "<<sq<<" not found!\n";
        abort();
      }
      break;
    case 'a':
      method = optarg;
      if (strcmp(method,"logh")&&strcmp(method,"ttl")) {
        std::cerr<<"Regularization method "<<method<<" not found!\n";
        abort();
      }
      break;
    case 'd':
      dtmin = atof(optarg);
      break;
    case 'e':
      err = atof(optarg);
      break;
    case 'p':
      pre = atoi(optarg);
      break;
    case 'h':
      std::cout<<"chain [option] filename\n"
               <<"Options: (*) show defaulted values\n"
               <<"    -N:  total number of particles ("<<n<<")\n"
               <<"    -n:  number of integration steps ("<<nstep<<")\n"
               <<"    -s:  step size, not physical time step ("<<s<<")\n"
               <<"    -a:  algorithmic regularization method; 'logh': Logarithmic Hamitonian; 'ttl': Time-transformed Leapfrog (logh)\n"
               <<"    -m:  use extrapolation method to get high accuracy; 'linear' for Romberg linear interpolation method; 'ration' for rational interpolation method (not used)\n"
               <<"    -M:  extrapolation sequences; 'even' for even sequences {h, h/2, h/4, h/8 ...}; 'bs' for Bulirsch & Stoer sequence {h, h/2, h/3, h/4, h/6, h/8 ...} (bs)\n"
               <<"    -e:  phase error limit ("<<err<<")\n"
               <<"    -d:  minimum physical time step ("<<dtmin<<")\n"
               <<"    -i:  maximum iteration steps for extrapolation method ("<<itermax<<")\n"
               <<"    -k:  sub-step number if no extrapolation method is used ("<<nsubstep<<")\n"
               <<"    -w:  print width of value ("<<w<<")\n"
               <<"    -p:  print digital precision ("<<pre<<")\n"
               <<"Input file format: each line: mass, x, y, z, vx, vy, vz\n";
      return 0;
    default:
      std::cerr<<"Unknown argument. check '-h' for help.\n";
      abort();
    }

  if (argc==1) {
    std::cerr<<"Please provide particle data filename\n";
    abort();
  }
  char* filename = argv[argc-1];      
  
  std::fstream fs;
  fs.open(filename,std::fstream::in);
  if(!fs.is_open()) {
    std::cerr<<"Error: Filename "<<filename<<" not found\n";
    abort();
  }

  chain<Particle> c(n);
  Particle *p=new Particle[n];
  for (int i=0;i<n;i++) {
    double x,y,z,vx,vy,vz,m;
    fs>>m>>x>>y>>z>>vx>>vy>>vz;
    p[i]=Particle(m,x,y,z,vx,vy,vz);
  }
  c.addP(n,p);
  
  c.set_abg(1.0,0.0,0.0);
  if (method)
    if (strcmp(method,"ttl")==0) c.set_abg(0.0,1.0,0.0,0.001);
  
  c.init(0.0);
  std::cout<<std::setprecision(pre);
  int* stepcount=new int[itermax+1];
  memset(stepcount,0,(itermax+1)*sizeof(int));
  
  //print
  std::cout<<"Time"
           <<std::setw(w)<<"E_err"
           <<std::setw(w)<<"Ekin"
           <<std::setw(w)<<"Pot"
           <<std::setw(w)<<"B"
           <<std::setw(w)<<"w"
           <<std::setw(w)<<"W"    
           <<std::setw(w)<<" "<<"mass-x-y-z-vx-vy-vz-for-each-particles"<<std::endl;
  for (int i=0;i<nstep;i++) {
    std::cout<<c.getTime()
             <<std::setw(w)<<c.getEkin()+c.getPot()+c.getB()
             <<std::setw(w)<<c.getEkin()
             <<std::setw(w)<<c.getPot()
             <<std::setw(w)<<c.getB()
             <<std::setw(w)<<c.getw()
             <<std::setw(w)<<c.getW();
    for (int j=0;j<n;j++) {
      std::cout<<std::setw(w)<<p[j].getMass();
      for (int k=0;k<3;k++) {
        std::cout<<std::setw(w)<<p[j].getPos()[k];
      }
      for (int k=0;k<3;k++) {
        std::cout<<std::setw(w)<<p[j].getVel()[k];
      }
    }
    std::cout<<std::endl;
    int icount=0;
    int msq=2;
    if (sq) 
      if (strcmp(sq,"even")==0) msq=1;
    if (sw==NULL) c.Leapfrog_step_forward(s,nsubstep,NULL,dtmin);
    else if (strcmp(sw,"linear")==0) icount = c.extrapolation_integration(s,err,itermax,1,msq,NULL,dtmin);
    else if (strcmp(sw,"rational")==0) icount = c.extrapolation_integration(s,err,itermax,2,msq,NULL,dtmin);
    stepcount[icount]++;
  }

  if (sw!=NULL) {
    std::cerr<<"Step histogram:\n I\tCount\n";
    for (int i=1;i<=itermax;i++) {
      std::cerr<<i<<"\t"<<stepcount[i]<<std::endl;
    }
  }

  return 0;
}
  
