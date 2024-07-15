include("GetMetric.jl");

function GetRayDirection(a,pos,fov,npix,ii,jj)
    function WrapAngle(a1,a2)
        x2=atan(a1,a2);
        x1=2.0*atan(0.5*sqrt(a1^2+a2^2));
        x1,x2
    end
    t,r,θ,ϕ=pos;
    
    gd=gdown(a,pos);
    g00=gd[1,1];
    g11 =gd[2,2];
    g22 =gd[3,3];
    g33=gd[4,4];
    g13=gd[1,4];
   
    e0=sqrt(-(g33/(g00*g33 - g13^2)))*[1.0;0.0;0.0;-(g13/g33)]';
    e1= [0.0;-(1/sqrt(g11));0.0;0.0]'; 
    e2=[0.0;0.0;1/sqrt(g22);0.0]'; 
    e3=[0.0,0.0,0.0,-(1/sqrt(g33))]'; 
   
    xscr=2.0*tan(fov/2)*(ii - 0.5*(npix+1.0))/npix;
    yscr=2.0*tan(fov/2)*(jj - 0.5*(npix+1.0))/npix;
    Tx,Px=WrapAngle(xscr,yscr);
    k=1.0;
    taudot=k*(-e0+cos(Tx).*e1+sin(Tx)*cos(Px).*e2+sin(Tx)*sin(Px).*e3);
 
    momentum=gd*taudot'
end
