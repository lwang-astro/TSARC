#include <iostream>
#include <fstream>
#include <getopt.h>
#include <iomanip>
#include <cmath>
#include "Newtonian_acceleration.h"

int main(int argc, char **argv){
   typedef Float Float3[3];
   bool iflag=false;
   int num=1;
   int WIDTH=22;
   int PRECISION=15;

   int copt;
   while ((copt = getopt(argc, argv, "in:w:p:h")) != -1)
       switch (copt) {
       case 'i':
           iflag=true;
           break;
       case 'n':
           num=atoi(optarg);
           break;
       case 'w':
           WIDTH=atoi(optarg);
           break;
       case 'p':
           PRECISION=atoi(optarg);
           break;
       case 'h':
           std::cout<<"keplerorbit [option] datafilename\n"
                    <<"    pariticle data (two lines): m, x(1:3), v(1:3)\n"
                    <<"    orbital data (one line):    m1, m2, xc(1:3), vc(1:3), semi, ecc, angle(1:3), ecc_anomaly, true_anomaly, mean_anomaly, period (last 3 parameters are not used)\n"
                    <<"Options: (*) show defaulted values\n"
                    <<"   -i:        read particle data, output kepler orbit data (one line)\n"
                    <<"   -n [int]:  number of pairs(1)\n"
                    <<"   -w [int]:  print width(22)\n"
                    <<"   -p [int]:  print precision(15)\n"
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

   std::cout<<std::setprecision(PRECISION);

   Float m1,m2,mc;
   Float3 x1,v1,x2,v2,xc,vc;
   Float ax,ecc,angle[3],ecc_anomaly,mean_anomaly,true_anomaly,period;

   if(iflag) {
       for(int i=0; i<num; i++) {
           fs>>m1>>x1[0]>>x1[1]>>x1[2]>>v1[0]>>v1[1]>>v1[2];
           fs>>m2>>x2[0]>>x2[1]>>x2[2]>>v2[0]>>v2[1]>>v2[2];
           mc = m1+m2;
           for (int j=0; j<3; j++) {
               xc[j] = (x1[j]*m1+x2[j]*m2)/mc;
               vc[j] = (v1[j]*m1+v2[j]*m2)/mc;
           }
           for (int j=0; j<3; j++) {
               x1[j] -= xc[j];
               x2[j] -= xc[j];
               v1[j] -= vc[j];
               v2[j] -= vc[j];
           }

           if (fs.eof()) {
               std::cerr<<"Error: data file reach end when reading pairs (current loaded pair number is "<<i+1<<"; required pair number "<<num<<std::endl;
               abort();
           }
           Float3 dx,dv;
           for (int j=0; j<3; j++) {
               dx[j] = x2[j]-x1[j];
               dv[j] = v2[j]-v1[j];
           }
           NTA::calc_kepler_orbit_par(ax,period,ecc,angle,true_anomaly,ecc_anomaly,mean_anomaly,mc,dx,dv);
           std::cout<<std::setw(WIDTH)<<m1
                    <<std::setw(WIDTH)<<m2
                    <<std::setw(WIDTH)<<xc[0]
                    <<std::setw(WIDTH)<<xc[1]
                    <<std::setw(WIDTH)<<xc[2]
                    <<std::setw(WIDTH)<<vc[0]
                    <<std::setw(WIDTH)<<vc[1]
                    <<std::setw(WIDTH)<<vc[2]
                    <<std::setw(WIDTH)<<ax
                    <<std::setw(WIDTH)<<ecc
                    <<std::setw(WIDTH)<<angle[0]
                    <<std::setw(WIDTH)<<angle[1]
                    <<std::setw(WIDTH)<<angle[2]
                    <<std::setw(WIDTH)<<ecc_anomaly
                    <<std::setw(WIDTH)<<true_anomaly
                    <<std::setw(WIDTH)<<mean_anomaly
                    <<std::setw(WIDTH)<<period
                    <<std::endl;
       }
   }
   else {
       for(int i=0; i<num; i++) {
           fs>>m1>>m2>>xc[0]>>xc[1]>>xc[2]>>vc[0]>>vc[1]>>vc[2]>>ax>>ecc>>angle[0]>>angle[1]>>angle[2]>>ecc_anomaly>>true_anomaly>>mean_anomaly>>period;
           if (fs.eof()) {
               std::cerr<<"Error: data file reach end when reading pairs (current loaded pair number is "<<i+1<<"; required pair number "<<num<<std::endl;
               abort();
           }
           NTA::kepler_orbit_generator(x1,x2,v1,v2,m1,m2,ax,ecc,angle,ecc_anomaly);
           for (int j=0; j<3; j++) {
               x1[j] += xc[j];
               x2[j] += xc[j];
               v1[j] += vc[j];
               v2[j] += vc[j];
           }
           std::cout<<std::setw(WIDTH)<<m1
                    <<std::setw(WIDTH)<<x1[0]
                    <<std::setw(WIDTH)<<x1[1]
                    <<std::setw(WIDTH)<<x1[2]
                    <<std::setw(WIDTH)<<v1[0]
                    <<std::setw(WIDTH)<<v1[1]
                    <<std::setw(WIDTH)<<v1[2]
                    <<std::endl
                    <<std::setw(WIDTH)<<m2
                    <<std::setw(WIDTH)<<x2[0]
                    <<std::setw(WIDTH)<<x2[1]
                    <<std::setw(WIDTH)<<x2[2]
                    <<std::setw(WIDTH)<<v2[0]
                    <<std::setw(WIDTH)<<v2[1]
                    <<std::setw(WIDTH)<<v2[2]
                    <<std::endl;
       }
   }

   return 0;
}
