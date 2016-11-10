#include "AR.h"
#include "particle.h"
#include <iostream>
#include <fstream>
//#include <unistd.h>
#include <getopt.h>
#include <string.h>
#include <stdlib.h>
#include <iomanip>
#include <cmath>

int main(int argc, char **argv){
  typedef double double3[3];
  std::size_t n=3; //particle number
  int w=18; //print width
  int pre=10; //print digital precision
  int nstep=1000; // total step size
  int nsubstep=128;  // sub-step number if direct LF method is used
  int itermax=20;  //iteration maximum number for extrapolation methods
  char* sw=NULL;        // if not 'none', use extrapolation method, 'linear' for Romberg method; 'rational'  for rational interpolation method
  char* sq=NULL;        // extrapolation sequence, 'rom' for {h,h/2,h/4,h/8...}; 'bs' for {h,h/2,h/3,h/4,h/6,h/8...}; '4k' for {h/2, h/6, h/10, h/14 ...}; 'hm' for {h, h/2, h/3, h/4 ...}
  char* method=NULL;   // regularization methods, 'logh': Logarithmic Hamitonian; 'ttl': Time-transformed Leapfrog\n (logh)
  double err=1e-10; // phase error requirement
  double s=0.5;    // step size
  double dtmin=5.4e-20; // mimimum physical time step
  double t=0.0;    // initial physical time
  double tend=-1; // ending physical time
  int copt;
  bool fflag=false; //external force flag
  double3* f=NULL;  // external force vectors
  
  //adjust step size
  bool dsA=false;   
  int nlev=8;    // adjust iteration level concentration for auto adjust step
  //  int nreduce=1;  // adjust reduction for smoothing
  //  double smooth=2.0;   // adjust smooth factor;


  static struct option long_options[] = {
    {"t-start", required_argument, 0, 0},
    {"t-end", required_argument, 0, 't'},
    {"nsub", required_argument, 0, 0},
    {"adjust-iter",required_argument, 0, 0},
    //    {"adjust-reduce",required_argument, 0, 0},
    //    {"adjust-smooth",required_argument, 0, 0},
    {"AR-method",required_argument, 0, 'r'},
    {"extra-method",required_argument, 0, 'm'},
    {"extra-seq",required_argument, 0, 'q'},
    {"iter-max",required_argument, 0, 0},
    {"error",required_argument, 0, 'e'},
    {"dtmin",required_argument, 0, 'd'},
    {"print-width",required_argument, 0, 0},
    {"print-precision",required_argument, 0, 0},
    {"help",no_argument, 0, 'h'},
    {0,0,0,0}
  };
  
  int option_index;
  while ((copt = getopt_long(argc, argv, "N:n:t:s:ar:m:q:i:e:d:fh", long_options, &option_index)) != -1)
    switch (copt) {
    case 0:
#ifdef DEBUG
      std::cerr<<"option "<<option_index<<" "<<long_options[option_index].name<<" "<<optarg<<std::endl;
#endif
      switch (option_index) {
      case 0:
        t = atof(optarg);
        break;
      case 2:
        nsubstep = atoi(optarg);
        break;
      case 3:
        nlev = atoi(optarg);
        if (nlev <=1) {
          std::cerr<<"Iteration expected level should be larger than 1\n";
          abort();
        }
        break;
        /*
      case 4:
        nreduce = atoi(optarg);
        break;
      case 5:
        smooth = atof(optarg);
        if (smooth<0) {
          std::cerr<<"Adjust step smooth factor cannot be negative\n";
          abort();
        }
        break;*/
      case 10:
        w = atof(optarg);
        break;
      case 11:
        pre = atoi(optarg);
        break;
      default:
        std::cerr<<"Unknown option. check '-h' for help.\n";
        abort();
      }
      break;
    case 'N':
      n = atoi(optarg);
      break;
    case 'n':
      nstep = atoi(optarg);
      break;
    case 't':
      tend = atof(optarg);
      break;
    case 's':
      s = atof(optarg);
      break;
    case 'a':
      dsA = true;
      break;
    case 'r':
      method = optarg;
      if (strcmp(method,"logh")&&strcmp(method,"ttl")) {
        std::cerr<<"Regularization method "<<method<<" not found!\n";
        abort();
      }
      break;
    case 'm':
      sw = optarg;
      if (strcmp(sw,"linear")&&strcmp(sw,"rational")&&strcmp(sw,"none")) {
        std::cerr<<"Extrapolation method "<<sw<<" not found!\n";
        abort();
      }
      break;
    case 'q':
      sq = optarg;
      if (strcmp(sq,"rom")&&strcmp(sq,"bs")&&strcmp(sq,"4k")&&strcmp(sq,"hm")) {
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
    case 'f':
      fflag = true;
      break;
    case 'h':
      std::cout<<"chain [option] data_filename\n"
               <<"Input data file format: each line: mass, x, y, z, vx, vy, vz\n"
               <<"Options: (*) show defaulted values\n"
               <<"    -N [int]:     total number of particles ("<<n<<")\n"
               <<"    -n [int]:     number of integration steps ("<<nstep<<")\n"
               <<"          --t-start [double]:  initial physical time ("<<t<<")\n"
               <<"    -t [double]:  ending physical time; if set, -n will be invalid (unset)\n"
               <<"          --t-end (same as -t)\n"
               <<"    -s [double]:  step size, not physical time step ("<<s<<")\n"
               <<"          --nsub [int]:       sub-step number if no extrapolation method is used ("<<nsubstep<<")\n"
               <<"    -a :          adjust step size automatically s *= pow(0.5, max(iter-reduce,0)/smooth) (off)\n"
               <<"          --adjust-iter   [int]:     iteration expected level for auto-adjust step size mode ("<<nlev<<")\n"
        //               <<"          --adjust-reduce [int]:     reduiteration expected level for auto-adjust step size mode ("<<nreduce<<")\n"
        //               <<"          --adjust-smooth [double]:  iteration expected level for auto-adjust step size mode ("<<smooth<<")\n"
               <<"    -r [string]:  algorithmic regularization method (logh)\n"
               <<"                  'logh': Logarithmic Hamitonian\n"
               <<"                  'ttl': Time-transformed Leapfrog\n"
               <<"          --AR-method (same as -r)\n"
               <<"    -m [string]:  use extrapolation method to get high accuracy (rational)\n"
               <<"                  'linear':   Romberg linear interpolation method;\n"
               <<"                  'rational': rational interpolation method;\n"
               <<"                  'none':     no extrapolation\n"
               <<"          --extra-method (same as -m)\n"
               <<"    -q [string]: extrapolation sequences (bs)\n"
               <<"                  'rom': Romberg sequence {h, h/2, h/4, h/8 ...};\n"
               <<"                  'bs':  Bulirsch & Stoer sequence {h, h/2, h/3, h/4, h/6, h/8 ...}\n"
               <<"                  '4k':  4k linear sequence {h/2, h/6, h/10, h/14 ...}\n"
               <<"                  'hm':  Harmonic sequence {h, h/2, h/3, h/4 ...}"
               <<"          --extra-seq (same as -q)\n"
               <<"    -i [int]: maximum iteration steps for extrapolation method ("<<itermax<<")\n"
               <<"          --iter-max (same as -i)\n"
               <<"    -e [double]:  phase and energy error limit ("<<err<<")\n"
               <<"          --error (same as -e)\n"
               <<"    -d [double]: [double]:  minimum physical time step ("<<dtmin<<")\n"
               <<"          --dtmin (same as -d)\n"
               <<"    -f :          use constant external force for each particles (read fx, fy, fz after reading each particle data in data file)\n"
               <<"          --print-width [int]:     print width of value ("<<w<<")\n"
               <<"          --print-precision [int]:     print digital precision ("<<pre<<")\n"
               <<"    -h :          print option information\n"
               <<"          --help (same as -h)\n";
      return 0;
    default:
      std::cerr<<"Unknown argument. check '-h' for help.\n";
      abort();
    }

  if (argc==1) {
    std::cerr<<"Please provide particle data filename\n";
    abort();
  }

#ifdef DEBUG
  std::cerr<<"Options:\n"
           <<"N: "<<n<<std::endl
           <<"steps: "<<nstep<<std::endl
           <<"t-start: "<<t<<std::endl
           <<"t-end: "<<tend<<std::endl
           <<"step size: "<<s<<std::endl
           <<"sub steps: "<<nsubstep<<std::endl
           <<"adjust: "<<dsA<<std::endl
           <<"adjust-iter: "<<nlev<<std::endl
    //           <<"adjust-reduce: "<<nreduce<<std::endl
    //           <<"adjust-smooth: "<<smooth<<std::endl
           <<"AR-method: "<<(method?method:"logH")<<std::endl
           <<"extra-method: "<<(sw?sw:"rational")<<std::endl
           <<"extra-seq: "<<(sq?sq:"bs")<<std::endl
           <<"itermax: "<<itermax<<std::endl
           <<"error: "<<err<<std::endl
           <<"dtmin: "<<dtmin<<std::endl
           <<"extra-force: "<<fflag<<std::endl
           <<"print width: "<<w<<std::endl
           <<"print precision: "<<pre<<std::endl;
#endif
  
  char* filename = argv[argc-1];      
  
  std::fstream fs;
  fs.open(filename,std::fstream::in);
  if(!fs.is_open()) {
    std::cerr<<"Error: Filename "<<filename<<" not found\n";
    abort();
  }

  chainpars pars;
  if (method)
    if (strcmp(method,"ttl")==0) {
      pars.setM(0.003,false);
      pars.setabg(0.0,1.0,0.0);
    }

  int msq=2;
  if (sq) {
    if (strcmp(sq,"rom")==0) msq=1;
    else if (strcmp(sq,"bs")==0) msq=2;
    else if (strcmp(sq,"4k")==0) msq=3;
    else msq=4;
  }
  int ms=2;
  if (sw) {
    if (strcmp(sw,"linear")==0) ms=1;
    else if (strcmp(sw,"none")==0) ms=0;
  }
  pars.setEXP(err,dtmin,1e-12,itermax,ms,msq,nlev);

  chain<Particle> c(n,pars);
  Particle *p=new Particle[n];
  if (fflag) f=new double3[n];
  for (std::size_t i=0;i<n;i++) {
    double x,y,z,vx,vy,vz,m;
    fs>>m>>x>>y>>z>>vx>>vy>>vz;
    if (fflag) fs>>f[i][0]>>f[i][1]>>f[i][2];
    if (fs.eof()) {
      std::cerr<<"Error: data file reach end when reading particles (current loaded particle number is "<<i<<"; required N = "<<n<<std::endl;
      abort();
    }
    p[i]=Particle(m,x,y,z,vx,vy,vz);
  }
  c.addP(n,p);

  c.init(t);
  std::cout<<std::setprecision(pre);

  //print
  std::cout<<"Time"
           <<std::setw(w)<<"E_err"
           <<std::setw(w)<<"Ekin"
           <<std::setw(w)<<"Pot"
           <<std::setw(w)<<"B"
           <<std::setw(w)<<"w"
           <<std::setw(w)<<"W"    
           <<std::setw(w)<<" "<<"mass-x-y-z-vx-vy-vz-for-each-particles"<<std::endl;
  int i=0;

  // step 
  double ds = s;

  // flag for output
  bool flag_out=true;
  
  while (true) {
    if (flag_out) {
      std::cout<<c.getTime()
               <<std::setw(w)<<(c.getEkin()+c.getPot()+c.getB())/c.getB()
               <<std::setw(w)<<c.getEkin()
               <<std::setw(w)<<c.getPot()
               <<std::setw(w)<<c.getB()
               <<std::setw(w)<<c.getw()
               <<std::setw(w)<<c.getW();
      for (std::size_t j=0;j<n;j++) {
        std::cout<<std::setw(w)<<p[j].getMass();
        for (int k=0;k<3;k++) {
          std::cout<<std::setw(w)<<p[j].getPos()[k];
        }
        for (int k=0;k<3;k++) {
          std::cout<<std::setw(w)<<p[j].getVel()[k];
        }
      }
      std::cout<<std::endl;
    }

    if ((tend<0&&i==nstep)||(tend>0&&std::abs(c.getTime()-tend)<1e-6)) break;
    i++;

    if (ms) {
      double dsf=c.extrapolation_integration(ds,tend,f);
      if (dsf<0) {
        ds *= -dsf;
        flag_out=false;
      }
      else if (dsA) {
        /*double factor=std::max(0.0,(std::abs(icount-nlev)-nreduce)/smooth);
        s *=std::pow(0.5,(icount>nlev?factor:-factor));
        s = std::min(s,0.9);*/
        ds = s*dsf;
        flag_out = true;
#ifdef DEBUG        
        std::cerr<<"S: "<<i<<" "<<ds<< std::endl;
#endif
      }
      else flag_out = true;
    }
      
    else c.Leapfrog_step_forward(s,nsubstep,f);
#ifdef TIME_PROFILE
    std::cerr<<"Time profile: Step: "<<i<<"  Accelaration+Potential(s): "<<c.profile.t_apw<<"  Update_link(s): "<<c.profile.t_uplink<<"  Leap-frog(s): "<<c.profile.t_lf<<"  Extrapolation(s): "<<c.profile.t_ep<<"  Perturbation(s): "<<c.profile.t_pext<<std::endl;
#endif
  }
  

  if (ms) {
    int* stepcount = c.profile.stepcount;
    std::cerr<<"Step histogram:\n I\tCount\tSteps\n";
    int subsum=0;
    int stepeven=2;
    int stepodd=3;
    int step=1;
    if (msq==3) step=2;
    int stepsum=0;
    int itersum=0;
    for (int i=1;i<=itermax;i++) {
      stepsum +=step;
      std::cerr<<i<<"\t"<<stepcount[i]<<"\t"<<step<<"\t"<<stepsum<<std::endl;
      subsum += stepcount[i]*stepsum;
      itersum += std::max(stepcount[i]*(i-1),0);
      if (msq==1) step *=2;
      else if (msq==2) {
        if (i%2) {
          step = stepeven;
          stepeven *=2;
        }
        else {
          step = stepodd;
          stepodd *=2;
        }
      }
      else step += 4;
    }
    std::cerr<<"Sum-of-steps: "<<subsum<<" iter-of-steps: "<<itersum<<std::endl;
  }

  delete[] p;
  if (fflag) delete[] f;
  
  return 0;
}

