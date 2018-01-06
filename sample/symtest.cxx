#include<iostream>
#include<cmath>
#include<iomanip>
#include "symplectic.h"

#define PRINT_WIDTH 20

int main(int argc, char** argv){
    int n = 1;
    if(argc==2) {
        n = atoi(argv[1])/2;
    }
    std::cout<<"Symplectic integrator order: "<<2*n<<std::endl;

    const int k= n>0?std::pow(3,n-1)+1:(n==-3?8:(n==-4?16:0));
    if(k==0) {
        std::cerr<<"Unknown order provided n = "<<2*n<<std::endl;
        abort();
    }
    double coff_cd[k][2];
    SYM::symcumck coff_order[k];
    
    SYM::symplectic_cofficients(coff_cd, coff_order, n, k);

    std::cout<<std::setprecision(15);
    std::cout<<std::setw(PRINT_WIDTH)<<"c_k"
             <<std::setw(PRINT_WIDTH)<<"d_k"
             <<std::setw(PRINT_WIDTH)<<"cum(c_k)"
             <<std::setw(PRINT_WIDTH)<<"k"
             <<std::endl;
    for (int i=0; i<k; i++) {
        std::cout<<std::setw(PRINT_WIDTH)<<coff_cd[i][0]
                 <<std::setw(PRINT_WIDTH)<<coff_cd[i][1]
                 <<std::setw(PRINT_WIDTH)<<coff_order[i].cck
                 <<std::setw(PRINT_WIDTH)<<coff_order[i].index
                 <<std::endl;
    }

    return 0;
}
