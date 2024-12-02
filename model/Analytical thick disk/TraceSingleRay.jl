include("GetRayDirection.jl");
include("ZamoProjection.jl");
include("ParallelTransport.jl");
include("PlasmaEquation.jl");
include("RayEqns.jl")



function TraceSingleRay(ν0,a,pos,fov,npix,ii,jj)

rh = 1.0+sqrt(1.0-a^2);


function event1(u,t,integrator)
       u[1];
end


function event!(out,u,t,integrator)
  out[1]=1.01*rh-u[2];
  out[2]=u[2]-500.0;
end

function affect!(integrator,idx)
  if idx == 1
    igt=terminate!(integrator);
  elseif idx == 2
    igt=terminate!(integrator);
  end
  igt
end

function raytra(du,u,p,τ)
    a,pt,pϕ=p;
    du[1:6]=RayEqns(a,u[1:4],[pt;u[5:6];pϕ])
end

function alleqn(du,u,p,τ)
    a,pt,pϕ,ν0=p;
    du[1:6]=RayEqns(a,u[1:4],[pt;u[5:6];pϕ]);
    du[7:22]=-ParallelTransport(a,u[1:4],[pt;u[5:6];pϕ],u[7:22])+PlasmaEquation(a,u[1:4],[pt;u[5:6];pϕ],u[7:22],ν0);
end

#Obtaining initial parameters
t,r,θ,ϕ=pos;
pt,pr,pθ,pϕ=GetRayDirection(a,pos,fov,npix,ii,jj);
u0 = [t,r, θ, ϕ, pr, pθ];
cache =(a,pt,pϕ);
prob = ODEProblem(raytra,u0,(0.0,3000.0),cache);
cb =VectorContinuousCallback(event!,affect!,nothing,2,save_positions=(true,true));

#Backward raytracing
sol = solve(prob,DP5(), abstol=1e-10, reltol=1e-9,callback = cb);
out=sol.u[end];

#Forward 
x0=[out[1:4];-out[5:6];zeros(Float64,16,1)];
para=(a,-pt,-pϕ,ν0);

pbe = ODEProblem(alleqn,x0,(0.0,3000.0),para);
ccb =ContinuousCallback(event1, igt->terminate!(igt), save_positions=(true,true));
solu= solve(pbe,DP5(),abstol=1e-7,reltol=1e-7,callback =ccb);

su=solu.u[end];

Λ=ZamoProjection(a,su[1:4]);

Sc=[su[7:10];;
[su[8];su[11:13]];;
[su[9];su[12];su[14:15]];;
[su[10],su[13],su[15],su[16]]
]';

sc=[[0.0,su[17],su[18],su[19]];;
[-su[17],0.0,su[20],su[21]];;
[-su[18],-su[20],0.0,su[22]];;
[-su[19],-su[21],-su[22],0.0]
]';

ted=Λ*Sc*(Λ');
trp=Λ*sc*(Λ');

II=0.5*(ted[2,2]+ted[3,3]);
QQ=0.5*(ted[2,2]-ted[3,3]);
UU=0.5*(ted[2,3]+ted[3,2]);
VV=0.5*(trp[2,3]-trp[3,2]);

[II,QQ,UU,VV]

end

