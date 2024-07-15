include("GetRadiationParameter.jl")
include("GetMetric.jl")
include("RedShiftFactor.jl")
include("FluidTetrads.jl")
include("FluidParameter.jl")

function PlasmaEquation(a,x,p,S,ν0)

N11,N12,N13,N14,N22,N23,N24,N33,N34,N44,n12,n13,n14,n23,n24,n34=S;


ne=NumberDensity(a,x);
Te=Temperature(a,x);
g=RedShiftFactor(a,x,p);
B_up=GetMagneticField(a,x);
B,θB=LocalB(a,x,p,B_up);


JI,JQ,JU,JV,aI,aQ,aU,aV,rQ,rU,rV=GetRadiationParameter(ne,Te,B,θB,ν0,g)


S1=-[[0.0,0.0,0.0,0.0];;
[0.0,aI+aQ,aU+rV,0.0];;
[0.0,aU-rV,aI-aQ,0.0];;
[0.0,0.0,0.0,0.0]]';

S3=[[0.0,0.0,0.0,0.0];;
[0.0,rQ,rU-aV,0.0];;
[0.0,rU+aV,-rQ,0.0];;
[0.0,0.0,0.0,0.0]]';

J=[[0.0,0.0,0.0,0.0];;
[0.0,JI+JQ,JU,0.0];;
[0.0,JU,JI-JQ,0.0];;
[0.0,0.0,0.0,0.0]]';

j=[[0.0,0.0,0.0,0.0];;
[0.0,0.0,JV,0.0];;
[0.0,-JV,0.0,0.0];;
[0.0,0.0,0.0,0.0]]';


Λ=FluidTetrads(a,x,p);
BLj=Λ'*j*Λ;
BLJ=Λ'*J*Λ;


BLS1=Λ'*S1*Λ;
BLS3=Λ'*S3*Λ;

gd=gdown(a,x);


SM=[[N11,N12,N13,N14];;
    [N12,N22,N23,N24];;
    [N13,N23,N33,N34];;
    [N14,N24,N34,N44]]';

sM=[[0.0,n12,n13,n14];;
    [-n12,0.0,n23,n24];;
    [-n13,-n23,0.0,n34];;
    [-n14,-n24,-n34,0.0]]';


dS=BLJ+0.5.*(BLS1*gd*SM+SM*gd*(BLS1')-BLS3*gd*sM+sM*gd*(BLS3'));

gds=BLj+0.5.*(BLS1*gd*sM+sM*gd*(BLS1')+BLS3*gd*SM-SM*gd*(BLS3'));

dS=0.5.*(dS+dS');
gds=0.5.*(gds-gds');

[dS[1,:];dS[2,2:4];dS[3,3:4];dS[4,4];gds[1,2:4];gds[2,3:4];gds[3,4]];
end


