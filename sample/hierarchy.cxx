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
  int w,pre;
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
    const Float* x[2];
    const Float* v[2];
    Float m[2];
    for (std::size_t i=0; i<2; i++) {
      x[i]=c[i]->getPos();
      v[i]=c[i]->getVel();
      m[i]=c[i]->getMass();
    }
        
    Float ax,per,ecc,angle[3],true_anomaly,ecc_anomaly,mean_anomaly; 
    Float dx[3] = {x[1][0]-x[0][0], x[1][1]-x[0][1], x[1][2]-x[0][2]};
    Float dv[3] = {v[1][0]-v[0][0], v[1][1]-v[0][1], v[1][2]-v[0][2]};
    Float mt = m[0]+m[1];
    
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

    Float xcm[3]={(x[0][0]*m[0]+x[1][0]*m[1])/mt, 
                  (x[0][1]*m[0]+x[1][1]*m[1])/mt, 
                  (x[0][2]*m[0]+x[1][2]*m[1])/mt};
    Float vcm[3]={(v[0][0]*m[0]+v[1][0]*m[1])/mt, 
                   (v[0][1]*m[0]+v[1][1]*m[1])/mt, 
                   (v[0][2]*m[0]+v[1][2]*m[1])/mt};

    return Particle(mt,xcm,vcm);
}

void chain_print(const ARC::chain<Particle> &c, const Float ds, const int w, const int pre) {
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
  typedef Float Float3[3];

  Float s=0.5;    // step size
  int n=1000; // total step size
  int m=0;    // method
  int k=4;    // symplectic integrator order or extrapolation method
  int itermax=20; // maximum iteration for extrapolation method
  int method=1; // regularization methods, 1: Logarithmic Hamitonian; 2: Time-transformed Leapfrog\n (logh)
  int imethod=2; // interpolation method (1. linear; 2. rational)
  bool pkepler=false; // whether print kepler parameters
  bool coffin=false; // whether only include inner binary potential for transformation function

  int copt;
  while ((copt = getopt(argc, argv, "n:s:hm:r:k:i:c:p")) != -1)
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
               <<"Example for hiarch_order & branch_id:\n"
               <<"     -------------------------------------------------\n"
               <<"       hiarch order          branch id                \n"
               <<"           0                      0                   \n"
               <<"                                 / \\                  \n"
               <<"           1                    0   1                 \n"
               <<"                               / \\ / \\                \n"
               <<"           2                  0  1 2  3               \n"
               <<"     -------------------------------------------------\n"
               <<"     PS: if 1-0 has no children, 1-1 still have 2, 3 as children's branch_id\n"
               <<"Options: (*) show defaulted values\n"
               <<"    -n [int]:     number of integration steps ("<<n<<")\n"
               <<"    -s [Float]:   step size, not physical time step ("<<s<<")\n"
               <<"    -m [int]:     integration methods: 0: symplectic, 1: extrapolation ("<<m<<")\n"
               <<"    -r [int]:     regularization methods, 1: Logarithmic Hamitonian; 2: Time-transformed Leapfrog ("<<method<<")\n"
               <<"    -k [int]:     if symplectic, k is order, if extrapolation, k is extrapolation sequence (1: Romberg sequence; 2. BS sequence; 3. 4k sequence; 4. Harmonic sequence) ("<<k<<")\n"
               <<"    -i [int]:     interpolation method (1. linear; 2. rational) ("<<imethod<<")\n"
               <<"    -c      :     set transformation function cofficient to only include inner most binary potential \n"
               <<"    -p      :     print Kepler orbital parameters\n";
      return 0;
    case 'm':
      m = atoi(optarg);
      break;
    case 'r':
      method = atoi(optarg);
      if(method!=1&&method!=2) {
          std::cerr<<"Error: regularization method should be either 1 or 2, provided: "<<method<<"\n";
          abort();
      }
      break;
    case 'k':
      k = atoi(optarg);
      break;
    case 'i':
      imethod = atoi(optarg);
      break;
    case 'c':
      coffin = true;
      break;
    case 'p':
      pkepler = true;
      break;
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

  if(method==1) pars.setabg(1.0,0.0,0.0);
  else if(method==2) pars.setabg(0.0,1.0,0.0);

  if(m==1) {
      pars.setIterSeq(itermax,k,10);
      pars.setIntp(imethod);
  }
  else pars.setSymOrder(k);
  //pars.setErr(err,dtmin,terr);
  //pars.setIterSeq(itermax,msq,intpmax);
  //pars.setIntp(ms);
  //pars.setAutoStep(dsA,std::max(std::min(itermax-3,5),1),std::min(std::max(itermax-5,3),itermax));
  //pars.setIterConst(iterfix);

  // print chain pars
  pars.print(std::cerr);

  int N;
  fs>>N;
  Particle p[N];
  
  Float p2min=-1.0;
  int idmin=0,ibmin=0;
  ptree<Particle, print_pars> plist;
  for(int i=0;i<N-1;i++) {
    int id,ib;
    Float m1,m2,ax,ecc,angle[3],ecc_anomaly;
    fs>>id>>ib>>m1>>m2>>ax>>ecc>>angle[0]>>angle[1]>>angle[2]>>ecc_anomaly;
    if (fs.eof()) {
      std::cerr<<"Error: data file reach end when reading pairs (current loaded pair number is "<<i<<"; required pair number "<<N-1<<std::endl;
      abort();
    }
    
    Float3 x1,x2,v1,v2;
    NTA::kepler_orbit_generator(x1,x2,v1,v2,m1,m2,ax,ecc,angle,ecc_anomaly);
    Float p2 = ax*ax*ax/(m1+m2);
    if (p2min<0||p2min>p2) {
        p2min = p2;
        idmin = id;
        ibmin = ib;
    }
    
    Particle a(m1,x1,v1);
    Particle b(m2,x2,v2);
    if (coffin) {
        a.setCoff(0);
        b.setCoff(0);
    }
    bool flag=plist.link(id,ib,a,b,pshift);
    if (!flag) {
      std::cerr<<"Error: particle id "<<id<<", ib "<<ib<<" are inconsistent with global particle tree structure, cannot created pairs!\n";
      abort();
    }
  }

  Particle* pa=plist.getP(idmin,ibmin,0);
  Particle* pb=plist.getP(idmin,ibmin,1);
  pa->setCoff(pa->getMass());
  pb->setCoff(pb->getMass());

  int count=plist.collect_and_store(p,N);
  if (count<0) {
    std::cerr<<"Error: particle number mismatched particle tree!\n";
    abort();
  }

//// for debugging
//  for (int i=0; i<N; i++) {
//    std::cout<<"m "<<p[i].getMass()<<" c "<<p[i].getCoff()<<" x "<<p[i].getPos()[0]<<std::endl;
//  }

  NTA::Newtonian_pars Int_pars;

  // new chain class
  ARC::chain<Particle> c(N);
  
  //c.link_int_par(Int_pars);
  c.linkP(N,p);

  c.init(0.0,pars,&Int_pars);

  // printing data
  print_pars pw;
  pw.w=18; //width
  pw.pre=10; //precision
  
  chain_print(c,s,pw.w,pw.pre);
  if (pkepler) plist.pair_process(0,0,kepler_print,pw);
  std::cout<<std::endl;

  // step size
  const Float ds = s;

  for(int i=0;i<n;i++) {
      if(m==1) {
          Float dsf=c.extrapolation_integration<Particle, ARC::Float3, NTA::Newtonian_pars>(ds,pars,-1,&Int_pars);
          if (dsf==0) {
              c.info->ErrMessage(std::cerr);
              abort();
          }
      }
      else{
          c.Symplectic_integration<Particle, ARC::Float3, NTA::Newtonian_pars>(ds, pars, NULL, &Int_pars);
      }

    chain_print(c,ds,pw.w,pw.pre);
    if (pkepler) plist.pair_process(0,0,kepler_print,pw);
    std::cout<<std::endl;
#ifdef ARC_PROFILE
    c.profile.print(std::cerr,i);
#endif
  }

#ifdef ARC_PROFILE
  if (m==1) {
    int* step=new int[itermax+1];
    EP::seq_BS(step,itermax+1);
    int* stepcount = c.profile.stepcount;
    std::cerr<<"Step histogram:\n I\tCount\tNstep\tStepsum\n";
    int stepsum=0;
    int itersum=0;
    for (int i=1;i<=itermax;i++) {
      stepsum += step[i-1];
      std::cerr<<i<<"\t"<<stepcount[i]<<"\t"<<step[i-1]<<"\t"<<stepsum<<std::endl;
      itersum += std::max(stepcount[i]*(i-1),0);
    }
    std::cerr<<"Sum-of-steps: "<<c.profile.itercount<<" Sum-of-iterations: "<<itersum<<std::endl;
    delete [] step;

  }
#endif
  
  return 0;
}
