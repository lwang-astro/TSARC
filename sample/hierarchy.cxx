#include "AR.h"
#include "particle.h"
#include "Newtonian_acceleration.h"
#include <iostream>
#include <fstream>
#include <getopt.h>
#include "ptree.h"

/// print Kepler orbit parameters

struct print_pars{
public:
  double w,pre;
};

Particle pshift(const Particle &a, const Particle &ref) {
  return Particle(a.getMass(),
                  a.getPos()[0]+ref.getPos()[0],
                  a.getPos()[1]+ref.getPos()[1],
                  a.getPos()[2]+ref.getPos()[2],
                  a.getVel()[0]+ref.getVel()[0],
                  a.getVel()[1]+ref.getVel()[1],
                  a.getVel()[2]+ref.getVel()[2]);
}

Particle kepler_print(const std::size_t id, const std::size_t ib, Particle* c[2], print_pars& ppars){
    const double* x[2];
    const double* v[2];
    double m[2];
    for (std::size_t i=0; i<2; i++) {
      x[i]=c[i]->getPos();
      v[i]=c[i]->getVel();
      m[i]=c[i]->getMass();
    }
        
    double ax,per,ecc,angle[3],true_anomaly,ecc_anomaly,mean_anomaly; 
    double dx[3] = {x[1][0]-x[0][0], x[1][1]-x[0][1], x[1][2]-x[0][2]};
    double dv[3] = {v[1][0]-v[0][0], v[1][1]-v[0][1], v[1][2]-v[0][2]};
    double mt = m[0]+m[1];
    
    NTA::calc_kepler_orbit_par(ax,per,ecc,angle,true_anomaly,ecc_anomaly,mean_anomaly,mt,dx,dv);
    std::cout<<std::setw(ppars.w)<<id
             <<std::setw(ppars.w)<<ib
             <<std::setw(ppars.w)<<ax
             <<std::setw(ppars.w)<<ecc
             <<std::setw(ppars.w)<<per
             <<std::setw(ppars.w)<<angle[0]
             <<std::setw(ppars.w)<<angle[1]
             <<std::setw(ppars.w)<<angle[2]
             <<std::setw(ppars.w)<<ecc_anomaly
             <<std::setw(ppars.w)<<true_anomaly
             <<std::setw(ppars.w)<<mean_anomaly;

    double xcm[3]={(x[0][0]*m[0]+x[1][0]*m[1])/mt, 
                  (x[0][1]*m[0]+x[1][1]*m[1])/mt, 
                  (x[0][2]*m[0]+x[1][2]*m[1])/mt};
    double vcm[3]={(v[0][0]*m[0]+v[1][0]*m[1])/mt, 
                   (v[0][1]*m[0]+v[1][1]*m[1])/mt, 
                   (v[0][2]*m[0]+v[1][2]*m[1])/mt};

    return Particle(mt,xcm,vcm);
}

void chain_print(const ARC::chain<Particle> &c, const double ds, const double w, const double pre) {
  // printing digital precision
  std::cout<<std::setprecision(pre);

  std::cout<<c.getTime()
           <<std::setw(w)<<(c.getEkin()+c.getPot()+c.getPt())/c.getPt()
           <<std::setw(w)<<c.getEkin()
           <<std::setw(w)<<c.getPot()
           <<std::setw(w)<<c.getPt()
           <<std::setw(w)<<c.getw()
           <<std::setw(w)<<c.getW();
  const int n = c.getN();
  for (int j=0;j<n;j++) {
    std::cout<<std::setw(w)<<c.getP(j).getMass();
    for (int k=0;k<3;k++) {
      std::cout<<std::setw(w)<<c.getP(j).getPos()[k];
    }
    for (int k=0;k<3;k++) {
      std::cout<<std::setw(w)<<c.getP(j).getVel()[k];
    }
  }
  std::cout<<std::setw(w)<<ds;
}


int main(int argc, char **argv){
  typedef double double3[3];

  double s=0.5;    // step size
  int n=1000; // total step size

  int copt;
  while ((copt = getopt(argc, argv, "n:s:h")) != -1)
    switch (copt) {
    case 'n':
      n = atoi(optarg);
      break;
    case 's':
      s = atof(optarg);
      break;
    case 'h':
      std::cout<<"chain [option] data_filename\n"
               <<"Input data file format: hiarch_order branch_id mass1,mass2,semi,ecc,angle[3],ecc_anomaly\n"
               <<"Options: (*) show defaulted values\n"
               <<"    -n [int]:     number of integration steps ("<<n<<")\n"
               <<"    -s [double]:  step size, not physical time step ("<<s<<")\n";
      return 0;
    default:
      std::cerr<<"Unknown argument. check '-h' for help.\n";
      abort();
    }
  
  if (argc==1) {
    std::cerr<<"Please provide particle data filename\n";
    abort();
  }

  // data file name
  char* filename = argv[argc-1];

  // open data file
  std::fstream fs;
  fs.open(filename,std::fstream::in);
  if(!fs.is_open()) {
    std::cerr<<"Error: Filename "<<filename<<" not found\n";
    abort();
  }

  // chain controller
  ARC::chainpars pars;

  pars.setA(NTA::Newtonian_AW,NTA::Newtonian_extAcc,NTA::Newtonian_kepler_period);

  //pars.setErr(err,dtmin,terr);
  //pars.setIterSeq(itermax,msq,intpmax);
  //pars.setIntp(ms);
  //pars.setAutoStep(dsA,std::max(std::min(itermax-3,5),1),std::min(std::max(itermax-5,3),itermax));
  //pars.setIterConst(iterfix);

  int N;
  fs>>N;
  Particle p[N];
  
  ptree<Particle, print_pars> plist;
  for(int i=0;i<N-1;i++) {
    int id,ib;
    double m1,m2,ax,ecc,angle[3],ecc_anomaly;
    fs>>id>>ib>>m1>>m2>>ax>>ecc>>angle[0]>>angle[1]>>angle[2]>>ecc_anomaly;
    if (fs.eof()) {
      std::cerr<<"Error: data file reach end when reading pairs (current loaded pair number is "<<i<<"; required pair number "<<N-1<<std::endl;
      abort();
    }
    
    double3 x1,x2,v1,v2;
    NTA::kepler_orbit_generator(x1,x2,v1,v2,m1,m2,ax,ecc,angle,ecc_anomaly);
    
    Particle a(m1,x1,v1);
    Particle b(m2,x2,v2);
    bool flag=plist.link(id,ib,a,b,pshift);
    if (!flag) {
      std::cerr<<"Error: particle id "<<id<<", ib "<<ib<<" are inconsistent with global particle tree structure, cannot created pairs!\n";
      abort();
    }
  }

  int count=plist.collect_and_store(p,N);
  if (count<0) {
    std::cerr<<"Error: particle number mismatched particle tree!\n";
    abort();
  }

//// for debugging
//  for (int i=0; i<N; i++) {
//    std::cout<<"m "<<p[i].getMass()<<" x "<<p[i].getPos()[0]<<std::endl;
//  }

  NTA::Newtonian_pars Int_pars;

  // new chain class
  ARC::chain<Particle> c(N);
  
  //c.link_int_par(Int_pars);
  c.addP(N,p);

  c.init(0.0,pars,&Int_pars);

  // printing data
  print_pars pw;
  pw.w=18; //width
  pw.pre=10; //precision
  
  chain_print(c,0,pw.w,pw.pre);
  plist.pair_process(0,0,kepler_print,pw);
  std::cout<<std::endl;

  // step size
  double ds = s;

  for(int i=0;i<n;i++) {
      double dsf=c.extrapolation_integration<Particle, ARC::double3, NTA::Newtonian_pars>(ds,pars,-1,&Int_pars);
    if (dsf==0) {
      c.info->ErrMessage(std::cerr);
      abort();
    }
    chain_print(c,ds,pw.w,pw.pre);
    plist.pair_process(0,0,kepler_print,pw);
    std::cout<<std::endl;
  }
  
  return 0;
}
