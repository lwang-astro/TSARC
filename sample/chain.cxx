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
  char* sw=NULL;        // if not 'none', use extrapolation method, 'linear' for Romberg method; 'rational'  for rational interpolation method
  char* sq=NULL;        // extrapolation sequence, 'rom' for {h,h/2,h/4,h/8...}; 'bs' for {h,h/2,h/3,h/4,h/6,h/8...}
  char* method=NULL;   // regularization methods, 'logh': Logarithmic Hamitonian; 'ttl': Time-transformed Leapfrog\n (logh)
  double err=1e-8; // phase error requirement
  double s=0.5;    // step size
  double dtmin=5.4e-20; // mimimum physical time step
  double t=0.0;    // initial physical time
  int copt;

  while ((copt = getopt(argc, argv, "N:n:t:s:a:k:m:M:i:e:d:w:p:h")) != -1)
    switch (copt) {
    case 'N':
      n = atoi(optarg);
      break;
    case 'n':
      nstep = atoi(optarg);
      break;
    case 't':
      t = atof(optarg);
      break;
    case 's':
      s = atof(optarg);
      break;
    case 'a':
      method = optarg;
      if (strcmp(method,"logh")&&strcmp(method,"ttl")) {
        std::cerr<<"Regularization method "<<method<<" not found!\n";
        abort();
      }
      break;
    case 'k':
      nsubstep = atoi(optarg);
      break;
    case 'm':
      sw = optarg;
      if (strcmp(sw,"linear")&&strcmp(sw,"rational")&&strcmp(sw,"none")) {
        std::cerr<<"Extrapolation method "<<sw<<" not found!\n";
        abort();
      }
      break;
    case 'M':
      sq = optarg;
      if (strcmp(sq,"rom")&&strcmp(sq,"bs")) {
        std::cerr<<"Extrapolation sequences "<<sq<<" not found!\n";
        abort();
      }
      break;
    case 'i':
      itermax = atoi(optarg);
      break;
    case 'e':
      err = atof(optarg);
      break;
    case 'd':
      dtmin = atof(optarg);
      break;
    case 'w':
      w = atof(optarg);
      break;
    case 'p':
      pre = atoi(optarg);
      break;
    case 'h':
      std::cout<<"chain [option] data_filename\n"
               <<"Input data file format: each line: mass, x, y, z, vx, vy, vz\n"
               <<"Options: (*) show defaulted values\n"
               <<"    -N [int]:     total number of particles ("<<n<<")\n"
               <<"    -n [int]      number of integration steps ("<<nstep<<")\n"
               <<"    -t [double]:  initial physical time ("<<t<<"\n"
               <<"    -s [double]:  step size, not physical time step ("<<s<<")\n"
               <<"    -a [string]:  algorithmic regularization method; 'logh': Logarithmic Hamitonian; 'ttl': Time-transformed Leapfrog (logh)\n"
               <<"    -k [int]:     sub-step number if no extrapolation method is used ("<<nsubstep<<")\n"
               <<"    -m [string]:  use extrapolation method to get high accuracy (rational)\n"
               <<"                  'linear':   Romberg linear interpolation method;\n"
               <<"                  'rational': rational interpolation method;\n"
               <<"                  'none':     no extrapolation\n"
               <<"    -M [string]:  extrapolation sequences (bs)\n"
               <<"                  'rom': Romberg sequences {h, h/2, h/4, h/8 ...};\n"
               <<"                  'bs':   Bulirsch & Stoer sequence {h, h/2, h/3, h/4, h/6, h/8 ...}\n"
               <<"    -i [int]:     maximum iteration steps for extrapolation method ("<<itermax<<")\n"
               <<"    -e [double]:  phase error limit ("<<err<<")\n"
               <<"    -d [double]:  minimum physical time step ("<<dtmin<<")\n"
               <<"    -w [int]:     print width of value ("<<w<<")\n"
               <<"    -p [int]:     print digital precision ("<<pre<<")\n";
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
    if (fs.eof()) {
      std::cerr<<"Error: data file reach end when reading particles (current loaded particle number is "<<i<<"; required N = "<<n<<std::endl;
      abort();
    }
    p[i]=Particle(m,x,y,z,vx,vy,vz);
  }
  c.addP(n,p);
  
  if (method)
    if (strcmp(method,"ttl")==0) c.setPars(0.0,1.0,0.0);
  
  c.init(t);
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
      if (strcmp(sq,"rom")==0) msq=1;
    if (sw==NULL) icount = c.extrapolation_integration(s,err,itermax,2,msq,NULL,dtmin);
    else if (strcmp(sw,"rational")==0) icount = c.extrapolation_integration(s,err,itermax,2,msq,NULL,dtmin);
    else if (strcmp(sw,"linear")==0) icount = c.extrapolation_integration(s,err,itermax,1,msq,NULL,dtmin);
    else if (strcmp(sw,"none")==0) c.Leapfrog_step_forward(s,nsubstep,NULL,dtmin);
    stepcount[icount]++;
#ifdef TIME_PROFILE
    std::cerr<<"Time profile: Step: "<<i<<"  Accelaration+Potential(s): "<<c.getTP_apw()<<"  Update_link(s): "<<c.getTP_uplink()<<"  Leap-frog(s): "<<c.getTP_lf()<<"  Extrapolation(s): "<<c.getTP_ep()<<"  Perturbation(s): "<<c.getTP_pext()<<std::endl;
#endif
  }

  if (sw==NULL||strcmp(sw,"none")!=0) {
    std::cerr<<"Step histogram:\n I\tCount\n";
    for (int i=1;i<=itermax;i++) {
      std::cerr<<i<<"\t"<<stepcount[i]<<std::endl;
    }
  }
  return 0;
}
  
