# ARC
Algorithmic regularization chain (ARC) implemented in C++

The major file in this code is AR.h, which includes the ARC class "chain", its parameter controller "chainpars" and particle list class "chainlist". The user can use these classes to generate their own code with ARC integrators.

There are other two header files:
extrapolation.h: The functions used for interpolation and extrapolation
particle.h: a example particle class can be used for template class chain.

In doc directory, html and latex format of Doxygen documents for all these three header files can be found with detailed description of the API.

In sample directory, a sample code is implemented which can be used to integrate few body systems with ARC methods. Use 'Make', the executable file "chain" will be generated. Using './chain -h', the user can get the idea how to use this sample code.
