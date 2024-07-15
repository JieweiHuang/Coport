include("GetMetric.jl");
include("FourVelocity.jl")
include("GetMagneticField.jl")



function Rscale(a,rr,th)
    E0=1;
    function R(r,theta)
    2*r^3+((-1)+E0^2)*r^4+2*a^2*((-1)+E0^2)*r^2*cos(
      theta)^2+a^4*((-1)+E0^2)*cos(theta)^4+2*r*((-1)*
      a^2*((-1)+E0^2)*cos(theta)^4+(a*E0+(-1)*a*((-1)+
      E0^2)^(1/2)*sin(theta)^2)^2)
    end
    rh=1+sqrt(1-a^2);
    (R(rh,th)/R(rr,th));
end

function NumberDensity(a,x)
    t,r,θ,ϕ=x;
    thj=pi/2;
    sigma=0.2;
    nrh=exp(-(sin(θ)-sin(thj))^2/sigma^2);
    scale=Rscale(a,r,θ);
    nrh*sqrt(scale)*1e6;
end
    
function Temperature(a,x)
    t,r,θ,ϕ=x;
    z=20.0;
    scale=Rscale(a,r,θ);
    tt=scale^((1.0+z)/(2.0+z)/3.0);
    tt*1e12;
end


function LocalB(a,x,p,b)
    gp=gup(a,x);
    gd=gdown(a,x);
    u=FourVelocity(a,x);
    k=gp*p;

    B0=u'*gd*b;
    k0=u'*p;

    uB=b+B0.*u;
    uk=k+k0.*u;

    b_dot_k=uB'*gd*uk;
    b_dot_b=uB'*gd*uB;
    k_dot_k=uk'*gd*uk;
    θB=acos(b_dot_k./(sqrt(k_dot_k)*sqrt(b_dot_b)));

    
    θB=max(min(θB,pi),0.0);

    B=sqrt(b_dot_b);

    B,θB
end


