include("GetMetric.jl");

function FluidTetrads(a,x,p,fvp,d)
    gd=gdown(a,x);
    gp=gup(a,x);
    
    k=gp*p;
    k1,k2,k3,k4=k;
    
    d=d./(sqrt(abs(d'*gd*d)));
    d1,d2,d3,d4=d;
    
    u=fvp;
    U1,U2,U3,U4=u;
    
    Beta=d'*gd*u;
    Omega=-p'*u;
    CC=d'*p/Omega-Beta;
    dnorm=d'*gd*d;
    NN=sqrt(abs(dnorm+Beta^2-CC^2));
    
   
    e0=u;
    e3=(k./Omega-u);
    e2=(d+Beta*u-CC*e3)./NN;
    e1=(sqrt(-det(gd))/Omega/NN).*gp*[d4*k3*U2+(-1)*d3*k4*U2+(-1)*d4*k2*U3+d2*k4*U3+d3*k2*
    U4+(-1)*d2*k3*U4;(-1)*d4*k3*U1+d3*k4*U1+d4*k1*U3+(-1)*
    d1*k4*U3+(-1)*d3*k1*U4+d1*k3*U4;d4*k2*U1+(-1)*d2*k4*
    U1+(-1)*d4*k1*U2+d1*k4*U2+d2*k1*U4+(-1)*d1*k2*U4;(-1)*
    d3*k2*U1+d2*k3*U1+d3*k1*U2+(-1)*d1*k3*U2+(-1)*d2*k1*
    U3+d1*k2*U3];
    
    [e0;;e1;;e2;;e3]';
end


