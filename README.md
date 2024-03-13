# VentricularFibrosis

This repository contains the code used to generate simulations for the paper "Re-entry in models of cardiac ventricular tissue with scar represented as a Gaussian random field"

All code is written in C, except for the utility used to generate and sample Gaussian random fields, which is written in Matlab.

Each directory of simulation code includes a set of source code that will compile to a single executable using the gnu C compiler using:

gcc -o<executable> *.c -I./ -lm

To run a simulation, the executable must be placed in a directory that includes a file called DiffusionCoefficient.txt, which is a plain text file containing floating point numbers on a 400 x 400 grid, where each number represents the diffusion coefficient at a particular grid point. These files can be produced by the utility file MakePatchyScar_isthmus.m. The directory must also contain a subdirectory called STFfiles, whch is where files containing snapshots of transmembrane voltage are written.
