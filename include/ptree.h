#pragma once

//template<class particle>
//particle pshift(const particle &a, const particle &ref) {
//  return particle(a.getMass(),
//                  a.getPos()[0]+ref.getPos()[0],
//                  a.getPos()[1]+ref.getPos()[1],
//                  a.getPos()[2]+ref.getPos()[2],
//                  a.getVel()[0]+ref.getVel()[0],
//                  a.getVel()[1]+ref.getVel()[1],
//                  a.getVel()[2]+ref.getVel()[2]);
//}

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

/// Class to make hierarchical N-body systems
/*!
  Depend on templete class particle and pair processing function parameters proc_params
*/
template<class particle, class proc_params>
class ptree{
private:
  void* p[2];
  bool lp[2]; ///indicator whether the leaf is ptree (true) or not
  bool collflag;

public:

  /// typedef of pair processing function
  /*!
    @param[in] id: current tree depth (root is 0)
    @param[in] ib: current leaf index (count from left to right, total size is \f$ 2^{id} \f$
    @param[in] c: two particle pointer to the two pair members
    @param[in,out] pars: proc_params type parameters used in the processing function
    \return: new particle generated by the function (e.g. center-of-mass particle)
   */
  typedef particle (*pair_proc_function) (const std::size_t id, const std::size_t ib, particle* c[], proc_params& pars);

  /// typedef of particle shifting function used during split
  /*!
    Split function create new branch with two particles, when the two particles are stored, they are applied ths shifting function by referring the old particle stored in the leaf
    @param[in] a: particle need to be shifted
    @param[in] ref: reference paraticle
    \return: new particle
   */
  typedef particle (*particle_shift) (const particle &a, const particle &ref);

  /// construct the tree py filling two particles 
  /*!
    @param[in] a: particle one to left leaf
    @param[in] b: particle two to right leaf
   */
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
      
  /// fill particles to two leaves/branches
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

  /// Split (delete) one leaf, create a new branch (ptree) and store two particles
  /*!
    @param[in] i: leaf index (0: left; 1: right; others return false)
    @param[in] a: left paricle
    @param[in] b: right particle
    @param[in] pshift: particle shifting function applied to the two particles referring to the origin particle before storing
    \return If the branch are successfully created and filled return true, otherwise return false
   */
    bool split(const std::size_t i, const particle &a, const particle &b, particle_shift pshift) {
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

  /// add a particle pair to one of the leaf
  /*!
    Add particle pair with depth #id and branch index #ib
    Example:
    -------------------------------------------------\\
      depth id              branch id                \\
          0                      0                   \\
                                / \                  \\
          1                    0   1                 \\
                              / \ / \                \\
          2                  0  1 2  3               \\
    -------------------------------------------------\\

    @param[in] id: depth of the tree, top is 0
    @param[in] ib: leaf index, counting from 0 from left to right (maximum index \f$ 2^{id} \f$) 
    @param[in] a: particle one
    @param[in] b: particle two
    @param[in] pshift: particle shifting function applied to the two particles referring to the origin particle stored at the splitted leaf. The two new particles generated will be stored
    \return true: successful adding
   */
    bool link(const std::size_t id, const std::size_t ib, const particle &a, const particle &b, particle_shift pshift) {
    if(id>1) {
      if(lp[ib/id]) return ((ptree<particle, proc_params>*)(this->p[ib/id]))->link(id-1,ib%id,a,b,pshift);
      else return false;
    }
    else if(id==1) return this->split(ib,a,b,pshift);
    else return this->fill(a,b);
  }

  /// collect particles into plist and map addresses to it
  /*! Scan tree and push back particles into plist, then delete the branch particle and link the corresponding particle address in plist. This means the data memory is moved from tree leafs to plist array.
    @param[in] plist: particle type data array
    @param[in] n: plist array size (maximum particle number that can be stored)
    \return remaining empty number in plist (if return value is -1: collection failed)
   */
  int collect_and_store(particle plist[], const int n) {
    if(!collflag) {
      collflag=true;
      if (n<0) return n;
      int k=n; //remaining number
      for (std::size_t i=0; i<2; i++) {
        if(lp[i]) k=((ptree<particle, proc_params>*)(this->p[i]))->collect_and_store(&plist[n-k],k);
        else if(k>0) {
          plist[n-k]=*(particle*)p[i];
          delete (particle*)p[i];
          p[i] = &plist[n-k];
          k--;
        }
        else return -1;
      }
      return k;
    }
    else return -1;
  }

  /// collect particles into plist
  /*! Scan tree and store particle address into plist. 
    @param[in] plist: particle address array
    @param[in] n: plist array size (maximum particle number that can be stored)
    \return remaining empty number in plist
   */
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

  /// Apply functions to pairs in the trees
  /*!
    
   */
  particle pair_process(const std::size_t id, const std::size_t ib, pair_proc_function f, proc_params &pars) {
    particle* c[2];
    bool newflag[2]={false};
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

  /// find particle
  /*! Get particle due to the index
    @param[in] id: current tree depth (root is 0)
    @param[in] ib: current leaf index (count from left to right, total size is \f$ 2^{id} \f$
    @param[in] im: leaf index 0: left, 1: right
   */
    particle *getP(const std::size_t id, const std::size_t ib, const std::size_t im) {
        if(id>0) {
            if(lp[ib/id]) return ((ptree<particle, proc_params>*)(this->p[ib/id]))->getP(id-1,ib%id,im);
            else return NULL;
        }
        else return (particle*)(this->p[im]);
    }

};


  
