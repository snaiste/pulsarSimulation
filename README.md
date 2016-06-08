### Pulsar Simulation

Command line compilation: g++ sim.cpp -o sim -lm -lfftw3 -std=c++0x

Requirements: [fftw3 library](http://www.fftw.org)

This simulation was developed to be used with the software that reads two’s-complement 8-bit (whole numbers in the range [-127; 128]), complex-sampled, dual circular polarisation, single-sub-band data from a specified file of DADA format. A file of DADA format contains 4096 byte ASCII header and complex data samples alternating like so: Real X, Imaginary X, Real Y, Imaginary Y, Real X, Imaginary X, etc. In this project, X and Y corresponded to right and left circular polarisation components respectively.

The simulation generates a complex-sampled, broadband, random white noise modulated by a Gaussian (or square) pulsar pulse profile. White noise was selected to fully represent a natural signal, where the amount of power per unit bandwidth is independent of frequency. Random white noise (x(t)) is real-valued, hence to transform it to analytical, z(t) = x(t) + i y(t), Hilbert transform is applied. 

More details about analytical signals and Hilbert transforms available [here](http://astronomy.swin.edu.au/research/theses/wvanstratenthesis.pdf) [Accessed: June 2016].

