# Coport
**Coport** is a Julia-based program designed for computing covariant polarized radiation transfer in any spacetime. This code is particularly useful for imaging black hole accretion systems. 

This document outlines the key components of **Coport** and provides detailed instructions for its use.

# Prerequisites
Before running the **Coport** code, ensure the following packages are installed:

1. `DifferentialEquations.jl`: For solving differential equations.
2. `LinearAlgebra.jl`: For performing linear algebra computations. 
3. `Interpolations.jl`: For interpolation calculations.
4. `MAT.jl`: For importing and exporting `.mat` files.

# Functions
The functions in the **Coport** program are named to reflect their specific functionalities.

1. `TraceRay.jl`: Handles the parallel computation of all rays.
2. `TraceSingleRay.jl`: Computes a single ray.
3. `GetMetric.jl`: Obtains the metric $g_{\mu \nu}$ and $g^{\mu \nu}$ .
4. `RayEqns.jl`: Contains the equations governing ray motion.
5. `GetRayDirection.jl`: Derives initial ray directions using a specific camera model.
6. `ParallelTransport.jl`: Manages the part of the equation used for parallel transport of the polarization tensor.
7. `PlasmaEquation.jl`: Describes emission, absorption, and Faraday rotation in the equations.
8. `FluidParameter.jl`: Defines various parameters of the fluid.
9. `GetRadiationParameter.jl`: Obtains covariant emission, absorption, and Faraday rotation coefficients.
10. `FluidTetrads.jl`: Obtains fluid tetrads.
11. `ZamoProjection.jl`: Projects the polarization tensor at the observer's screen.
12. `DataInterpolation.jl`: Interpolates GRMHD data.

# Running
The settings for the observer and observation frequency in the **Coport** code are both located in the `main.jl` file. After configuring all parameters, you can simply run the `main.jl` file to generate the output. 

For instance, if we set the pixel density to $n\times n$, then a file named $\rm{Stokes}$ with dimensions $n\times n\times 4$ will be output. The Stokes parameters in this file will be $I=\rm{Stokes}[:,:,1],Q=\rm{Stokes}[:,:,2],U=\rm{Stokes}[:,:,3],V=\rm{Stokes}[:,:,4]$ , respectively.

# Credit
If you use our code, please cite
'[Coport: A New Public Code for Polarized Radiative Transfer in a Covariant Framework](https://arxiv.org/abs/2407.10431)'
by Jiewei Huang, Liheng Zheng, Minyong Guo and Bin Chen.
