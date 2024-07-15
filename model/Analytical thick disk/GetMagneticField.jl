include("FourVelocity.jl");

function GetMagneticField(a,x)
    t,r,θ,ϕ=x;
    MetricDown=gdown(a,x);
    g=sqrt(-det(MetricDown));
    fvp=FourVelocity(a,x);
    fvd=MetricDown*fvp;
    F=-sign(cos(θ))*sin(θ)*(10.0/0.447742572591883);
    output=F/(g*fvp[2]).*([1,0,0,0]+fvd[1].*fvp);
end
