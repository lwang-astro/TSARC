# ARC
Algorithmic regularization chain (ARC) implemented in C++

The major file in this code is AR.h, which includes the ARC class "chain", its parameter controller "chainpars" and particle list class "chainlist". The user can use these classes to generate their own code with ARC integrators.

There are other two header files:
extrapolation.h: The functions used for interpolation and extrapolation
particle.h: a example particle class can be used for template class chain.

In doc directory, html and latex format of Doxygen documents for all these three header files can be found with detailed description of the API.

In sample directory, two sample codes are provided.
The 'chain' program can read particle data and integrate the orbits using ARC library.
The 'hierarchy' program can read kepler orbital parameters and generate hierarchical few-body systems and do integration using ARC.
Using 'Make chain' and 'Make hierarchy', these programs can be compiled individually.
For the details of how to use them, using './chain -h' or './hierarchy -h'.
