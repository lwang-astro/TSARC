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

    const int k= std::pow(3,n-1)+1;
    double coff_cd[k][2];
    
    SYM::symplectic_cofficients(coff_cd, n, k);

    std::cout<<std::setprecision(15);
    std::cout<<std::setw(PRINT_WIDTH)<<"c_k"
             <<std::setw(PRINT_WIDTH)<<"d_k"
             <<std::endl;
    for (int i=0; i<k; i++) {
        std::cout<<std::setw(PRINT_WIDTH)<<coff_cd[i][0]
                 <<std::setw(PRINT_WIDTH)<<coff_cd[i][1]
                 <<std::endl;
    }

    return 0;
}
