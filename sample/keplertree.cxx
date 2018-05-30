#include <iostream>
#include <fstream>
#include <getopt.h>
#include <iomanip>
#include <cmath>
#include <vector>
#include "ptree.h"
#include "Newtonian_acceleration.h"

Particle pshift(const Particle &a, const Particle &ref) {
  return Particle(a.getMass(),
                  a.getPos()[0]+ref.getPos()[0],
                  a.getPos()[1]+ref.getPos()[1],
                  a.getPos()[2]+ref.getPos()[2],
                  a.getVel()[0]+ref.getVel()[0],
                  a.getVel()[1]+ref.getVel()[1],
                  a.getVel()[2]+ref.getVel()[2]);
}

int main(int argc, char **argv){
    typedef Float Float3[3];
    int num=0;
    int WIDTH=22;
    int PRECISION=15;
    int unit=0;

    int copt;
    while ((copt = getopt(argc, argv, "n:w:p:u:h")) != -1)
        switch (copt) {
        case 'n':
            num=atoi(optarg);
            break;
        case 'w':
            WIDTH=atoi(optarg);
            break;
        case 'p':
            PRECISION=atoi(optarg);
            break;
        case 'u':
            unit=atoi(optarg);
            break;
        case 'h':
            std::cout<<"keplertree [option] datafilename\n"
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
                     <<"   -n [int]:  number of pairs to read (defaulted: all)\n"
                     <<"   -w [int]:  print width(22)\n"
                     <<"   -p [int]:  print precision(15)\n"
                     <<"   -u [int]:  0: unscale; 1: x[PC], v[km/s], semi[AU], period[days]; 2: x[AU], v[km/s], semi[AU], period[days]\n"
                     <<std::endl;
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

    ptree<Particle, int> plist;
    int N=0;
    while (!fs.eof()) {
        int id,ib;
        Float m1,m2,ax,ecc,angle[3],ecc_anomaly;
        fs>>id>>ib>>m1>>m2>>ax>>ecc>>angle[0]>>angle[1]>>angle[2]>>ecc_anomaly;
        if (fs.eof()) break;
        N++;

        Float3 x1,x2,v1,v2;
        NTA::kepler_orbit_generator(x1,x2,v1,v2,m1,m2,ax,ecc,angle,ecc_anomaly,unit);
    
        Particle a(m1,x1,v1);
        Particle b(m2,x2,v2);
        bool flag=plist.link(id,ib,a,b,pshift);
        if (!flag) {
            std::cerr<<"Error: particle id "<<id<<", ib "<<ib<<" are inconsistent with global particle tree structure, cannot created pairs!\n";
            abort();
        }
        if (N==num) break;
    }

    N++;
    Particle p[N];

    int count=plist.collect_and_store(p,N);
    if (count<0) {
        std::cerr<<"Error: particle number mismatched particle tree!\n";
        abort();
    }

    std::cout<<std::setprecision(PRECISION);
    for (int i=0; i<N; i++) {
        std::cout<<std::setw(WIDTH)<<p[i].getMass()
                 <<std::setw(WIDTH)<<p[i].getPos()[0]
                 <<std::setw(WIDTH)<<p[i].getPos()[1]
                 <<std::setw(WIDTH)<<p[i].getPos()[2]
                 <<std::setw(WIDTH)<<p[i].getVel()[0]
                 <<std::setw(WIDTH)<<p[i].getVel()[1]
                 <<std::setw(WIDTH)<<p[i].getVel()[2]
                 <<std::endl;
    }
    
    return 0;
}
