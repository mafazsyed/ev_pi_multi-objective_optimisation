# Optimized PI Controller Tuning forElectric Vehicle using NSGA-II
Contains implementation of multi-objective optisation for PI turning in EVs using NSGA-II algorithm.

This repository contains components and metrics essential for computational performance, developed in C and accessible in MATLAB via wrappers. To use these functions in MATLAB, they must be compiled for your specific system.

## Compilation Instructions

### Compile Standard Functions
For standard functions like the simulated binary crossover operator, use the following MATLAB command:
```matlab
mex sbx.c

### Compile Hypervolume Metric
The hypervolume metric compilation is more involved and requires the following command:
mex -DVARIANT=4 Hypervolume MEX.c hv.c avl.c
