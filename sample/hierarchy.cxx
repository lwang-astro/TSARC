#include "AR.h"
#include "particle.h"
#include "Newtonian_acceleration.h"
#include <iostream>
#include <fstream>
#include <getopt.h>

Particle pshift(const Particle &a, const Particle &ref) {
  return Particle(a.getMass(),
                  a.getPos()[0]+ref.getPos()[0],
                  a.getPos()[1]+ref.getPos()[1],
                  a.getPos()[2]+ref.getPos()[2],
                  a.getVel()[0]+ref.getVel()[0],
                  a.getVel()[1]+ref.getVel()[1],
                  a.getVel()[2]+ref.getVel()[2]);
}

//Particle center_of_mass(const Particle &a, const Particle &b) {
//  double m1 = a.getMass();
//  double m2 = b.getMass();
//  double m = m1+m2;
//  double* x1 = a.getPos();
//  double* x2 = b.getPos();
//  double* v1 = a.getVel();
//  double* v2 = b.getVel();
//  return Particle(m,
//                  (x1[0]*m1+x2[0]*m2)/m,
//                  (x1[1]*m1+x2[1]*m2)/m,
//                  (x1[2]*m1+x2[2]*m2)/m,
//                  (v1[0]*m1+v2[0]*m2)/m,
//                  (v1[1]*m1+v2[1]*m2)/m,
//                  (v1[2]*m1+v2[2]*m2)/m);
//}
            

struct ptree{
private:
  void* p[2];
  bool lp[2];
  bool collflag;

  
public:
  ptree(const Particle& a, const Particle& b): collflag(false) {
    fill(a,b);
  }

  ptree(): collflag(false) {p[0]=p[1]=NULL; lp[0]=lp[1]=false;}

  ~ptree() {clear();}

  void clear() {
    for (std::size_t i=0; i<2; i++) {
      if(lp[i]) {
        ((ptree*)p[i])->clear();
        lp[i]=false;
      }
      else if(p[i]!=NULL) {
        if(!collflag) delete (Particle*)p[i];
        p[i]=NULL;
      }
    }
  }
      
  /// fill particles to two branches
  /*!
    @param[in] a: left particle
    @param[in] b: right particle
    \return If the branches are already filled (failure), return false, else true
   */
  bool fill(const Particle &a, const Particle &b) {
    if(p[0]==NULL) p[0]=new Particle(a);
    else return false;
    if(p[1]==NULL) p[1]=new Particle(b);
    else return false;
    lp[0]=lp[1]=false;
    return true;
  }

  bool split(const std::size_t i, const Particle &a, const Particle &b) {
    if(i<0||i>1) return false;
    if(!lp[i]) {
      if(p[i]!=NULL) {
        Particle* tmp=(Particle*)p[i];
        p[i]=new ptree;
        bool fg=((ptree*)p[i])->fill(pshift(a,*tmp),pshift(b,*tmp));
        delete tmp;
        if (!fg) return false;
      }
      else {
        p[i]=new ptree;
        bool fg=((ptree*)p[i])->fill(a,b);
        if (!fg) return false;
      }
      lp[i]=true;
      return true;
    }
    else return false;
  }

  /// add particle pair at the ending branch
  /*!
    @param[in] id: level of the tree, top level is 0
    @param[in] ib: branch index, counting from 0 from left to right
    @param[in] a: particle one
    @param[in] b: particle two
    \return true: successful adding
   */
  bool link(const std::size_t id, const std::size_t ib, const Particle &a, const Particle &b) {
    if(id>1) {
      if(lp[ib/id]) return ((ptree*)(this->p[ib/id]))->link(id-1,ib%id,a,b);
      else return false;
    }
    else if(id==1) return this->split(ib,a,b);
    else return this->fill(a,b);
  }

  /// collect particles into plist and link to it
  /*! Scan tree and push back particles into plist, then delete the branch particle and link the corresponding particle address in plist. Thus this function is used to moving particle data memory from tree branch to plist array.
    @param[in] plist: particle type data array
    @param[in] n: plist array size (maximum particle number that can be stored)
    \return remaining empty number in plist
   */
  int collect(Particle plist[], const int n) {
    collflag=true;
    if (n<0) return n;
    int k=n; //remaining number
    for (std::size_t i=0; i<2; i++) {
      if(lp[i]) k=((ptree*)(this->p[i]))->collect(&plist[n-k],k);
      else if(k>0) {
        plist[n-k].set(*(Particle*)p[i]);
        delete (Particle*)p[i];
        p[i] = &plist[n-k];
        k--;
      }
      else return -1;
    }
    return k;
  }

  /// print Kepler orbit parameters
  /*!
    \return return the center-of-mass of the current tree pair
   */
  Particle kepler_print(const int id, const int ib, const double w, const double pre){
    const double* x[2];
    const double* v[2];
    Particle* c[2];
    bool newflag[2]={false};
    double m[2];
    for (std::size_t i=0; i<2; i++) {
      if(lp[i]) {
        c[i]=new Particle(((ptree*)p[i])->kepler_print(id+1, (i==0?ib:ib+id+1),w,pre));
        newflag[i]=true;
      }
      else c[i]=(Particle*)p[i];
      x[i]=c[i]->getPos();
      v[i]=c[i]->getVel();
      m[i]=c[i]->getMass();
    }
        
    double ax,per,ecc,angle[3],true_anomaly,ecc_anomaly,mean_anomaly; 
    double dx[3] = {x[1][0]-x[0][0], x[1][1]-x[0][1], x[1][2]-x[0][2]};
    double dv[3] = {v[1][0]-v[0][0], v[1][1]-v[0][1], v[1][2]-v[0][2]};
    double mt = m[0]+m[1];
    
    NTA::calc_kepler_orbit_par(ax,per,ecc,angle,true_anomaly,ecc_anomaly,mean_anomaly,mt,dx,dv);
    std::cout<<std::setw(w)<<id
             <<std::setw(w)<<ib
             <<std::setw(w)<<ax
             <<std::setw(w)<<ecc
             <<std::setw(w)<<per
             <<std::setw(w)<<angle[0]
             <<std::setw(w)<<angle[1]
             <<std::setw(w)<<angle[2]
             <<std::setw(w)<<ecc_anomaly
             <<std::setw(w)<<true_anomaly
             <<std::setw(w)<<mean_anomaly;

    double xcm[3]={(x[0][0]*m[0]+x[1][0]*m[1])/mt, 
                  (x[0][1]*m[0]+x[1][1]*m[1])/mt, 
                  (x[0][2]*m[0]+x[1][2]*m[1])/mt};
    double vcm[3]={(v[0][0]*m[0]+v[1][0]*m[1])/mt, 
                   (v[0][1]*m[0]+v[1][1]*m[1])/mt, 
                   (v[0][2]*m[0]+v[1][2]*m[1])/mt};
    for (int i=0;i<2;i++) if(newflag[i]) delete c[i];

    return Particle(mt,xcm,vcm);
  }

};

void chain_print(const ARC::chain<Particle,NTA::Newtonian_pars> &c, const double ds, const double w, const double pre) {
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
  ARC::chainpars<NTA::Newtonian_pars> pars;

  pars.setA(NTA::Newtonian_AW,NTA::Newtonian_Ap,NTA::Newtonian_kepler_period);

  //pars.setErr(err,dtmin,terr);
  //pars.setIterSeq(itermax,msq,intpmax);
  //pars.setIntp(ms);
  //pars.setAutoStep(dsA,std::max(std::min(itermax-3,5),1),std::min(std::max(itermax-5,3),itermax));
  //pars.setIterConst(iterfix);

  int N;
  fs>>N;
  Particle p[N];
  
  ptree plist;
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
    bool flag=plist.link(id,ib,a,b);
    if (!flag) {
      std::cerr<<"Error: particle id "<<id<<", ib "<<ib<<" are inconsistent with global particle tree structure, cannot created pairs!\n";
      abort();
    }
  }

  int count=plist.collect(p,N);
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
  ARC::chain<Particle,NTA::Newtonian_pars> c(N,pars);
  
  c.link_int_par(Int_pars);
  c.addP(N,p);

  c.init(0.0);

  // printing data
  chain_print(c,0,18,10);
  plist.kepler_print(0,0,18,10);
  std::cout<<std::endl;

  // step size
  double ds = s;

  for(int i=0;i<n;i++) {
    double dsf=c.extrapolation_integration(ds);
    if (dsf==0) {
      c.info->ErrMessage(std::cerr);
      abort();
    }
    chain_print(c,ds,18,10);
    plist.kepler_print(0,0,18,10);
    std::cout<<std::endl;
  }
  
  return 0;
}
