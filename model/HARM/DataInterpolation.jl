function GetGRMHD()
    
Rin=matread("data\\r.mat")["Rin"];
THin=matread("data\\theta.mat")["THin"];
Pin=matread("data\\phi.mat")["Pin"];
r=Rin[:,1,1];
θ=THin[1,1,:];
ϕ=Pin[1,:,1];

n=matread("data\\n.mat")["nn"];
p=matread("data\\p.mat")["p1"];

Tp=matread("data\\Tp.mat")["Tp"];


B1=matread("data\\B1.mat")["zB1"];
B2=matread("data\\B2.mat")["zB2"];
B3=matread("data\\B3.mat")["zB3"];

U0=matread("data\\U0.mat")["u0all"];
U1=matread("data\\U1.mat")["u1all"];
U2=matread("data\\U2.mat")["u2all"];
U3=matread("data\\U3.mat")["u3all"];

intn=linear_interpolation((r,ϕ,θ),n,extrapolation_bc= Line());
intp=linear_interpolation((r,ϕ,θ),p,extrapolation_bc= Line());

intTp=linear_interpolation((r,ϕ,θ),Tp,extrapolation_bc= Line());

intB1=linear_interpolation((r,ϕ,θ),B1,extrapolation_bc= Line());
intB2=linear_interpolation((r,ϕ,θ),B2,extrapolation_bc= Line());
intB3=linear_interpolation((r,ϕ,θ),B3,extrapolation_bc= Line());

intU0=linear_interpolation((r,ϕ,θ),U0,extrapolation_bc= Line());
intU1=linear_interpolation((r,ϕ,θ),U1,extrapolation_bc= Line());
intU2=linear_interpolation((r,ϕ,θ),U2,extrapolation_bc= Line());
intU3=linear_interpolation((r,ϕ,θ),U3,extrapolation_bc= Line());


intn,intp,intTp,intB1,intB2,intB3,intU0,intU1,intU2,intU3


end