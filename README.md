# VentricularFibrosis

This repository contains the code used to generate simulations for the paper "Re-entry in models of cardiac ventricular tissue with scar represented as a Gaussian random field"

All code is written in C, except for the utility used to generate and sample Gaussian random fields, which is written in Matlab.

Each directory of simulation code includes a set of source code that will compile to a single executable using the gnu C compiler using

gcc -o<executable> *.c -I./ -lm


