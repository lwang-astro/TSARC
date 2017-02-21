#pragma once

template<class particle>
particle pshift(const particle &a, const particle &ref) {
  return particle(a.getMass(),
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
            
template<class particle, class proc_params>
class ptree{
private:
  void* p[2];
  bool lp[2];
  bool collflag;

  
public:

  typedef particle (*pair_proc_function) (const std::size_t, const std::size_t, particle* c[], proc_params& pars);

  ptree(const particle& a, const particle& b): collflag(false) {
    fill(a,b);
  }

  ptree(): collflag(false) {p[0]=p[1]=NULL; lp[0]=lp[1]=false;}

  ~ptree() {clear();}

  void clear() {
    for (std::size_t i=0; i<2; i++) {
      if(lp[i]) {
        ((ptree<particle, proc_params>*)p[i])->clear();
        lp[i]=false;
      }
      else if(p[i]!=NULL) {
        if(!collflag) delete (particle*)p[i];
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
  bool fill(const particle &a, const particle &b) {
    if(p[0]==NULL) p[0]=new particle(a);
    else return false;
    if(p[1]==NULL) p[1]=new particle(b);
    else return false;
    lp[0]=lp[1]=false;
    return true;
  }

  bool split(const std::size_t i, const particle &a, const particle &b) {
    if(i<0||i>1) return false;
    if(!lp[i]) {
      if(p[i]!=NULL) {
        particle* tmp=(particle*)p[i];
        p[i]=new ptree<particle, proc_params>;
        bool fg=((ptree<particle, proc_params>*)p[i])->fill(pshift(a,*tmp),pshift(b,*tmp));
        delete tmp;
        if (!fg) return false;
      }
      else {
        p[i]=new ptree<particle, proc_params>;
        bool fg=((ptree<particle, proc_params>*)p[i])->fill(a,b);
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
  bool link(const std::size_t id, const std::size_t ib, const particle &a, const particle &b) {
    if(id>1) {
      if(lp[ib/id]) return ((ptree<particle, proc_params>*)(this->p[ib/id]))->link(id-1,ib%id,a,b);
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
  int collect(particle plist[], const int n) {
    collflag=true;
    if (n<0) return n;
    int k=n; //remaining number
    for (std::size_t i=0; i<2; i++) {
      if(lp[i]) k=((ptree<particle, proc_params>*)(this->p[i]))->collect(&plist[n-k],k);
      else if(k>0) {
        plist[n-k].set(*(particle*)p[i]);
        delete (particle*)p[i];
        p[i] = &plist[n-k];
        k--;
      }
      else return -1;
    }
    return k;
  }

  
  int collect(particle* plist[], const int n) {
    if (n<0) return n;
    int k=n; //remaining number
    for (std::size_t i=0; i<2; i++) {
      if(lp[i]) k=((ptree<particle, proc_params>*)(this->p[i]))->collect(&plist[n-k],k);
      else if(k>0) {
        plist[n-k]=(particle*)p[i];
        k--;
      }
      else return -1;
    }
    return k;
  }

  particle pair_process(const std::size_t id, const std::size_t ib, pair_proc_function f, proc_params &pars) {
    particle* c[2];
    bool newflag[2]={false};
    double m[2];
    for (std::size_t i=0; i<2; i++) {
      if(lp[i]) {
        c[i]=new particle(((ptree<particle, proc_params>*)p[i])->pair_process(id+1, (i==0?ib:ib+id+1), f,  pars));
        newflag[i]=true;
      }
      else c[i]=(particle*)p[i];
    }
    particle nc(f(id, ib, c, pars));
    
    for (int i=0;i<2;i++) if(newflag[i]) delete c[i];
    return nc;
  }

};


  
