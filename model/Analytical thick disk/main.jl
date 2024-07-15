using DifferentialEquations
using LinearAlgebra
using MAT
include("TraceRay.jl")


a=0.94;
pos=[0.0,600.0,17*pi/180,0.0]
fov=pi/96;
npix=128::Int64;
ν0=230e9;
@time out=TraceRay(ν0,a,pos,fov,npix);


#Export
file = matopen("Stokes.mat", "w")
write(file, "Stokes", out)
close(file)

