using MAT
using Interpolations
using DifferentialEquations
using LinearAlgebra

include("TraceRay.jl")
include("DataInterpolation.jl")

struct Foo
    n
    p
    Tp
    B1
    B2
    B3
    U0
    U1
    U2
    U3
end

intn,intp,intTp,intB1,intB2,intB3,intU0,intU1,intU2,intU3=GetGRMHD();

str=Foo(intn,intp,intTp,intB1,intB2,intB3,intU0,intU1,intU2,intU3);




a=0.94;
pos=[0.0,600.0,17*pi/180,0.0]
fov=pi/96;
npix=128::Int64;
ν0=230e9;
@time out=TraceRay(ν0,a,pos,fov,npix,str);



file = matopen("Stokes.mat", "w")
write(file, "Stokes", out)
close(file)
