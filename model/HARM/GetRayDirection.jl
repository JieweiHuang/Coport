include("GetMetric.jl");

function GetRayDirection(a,pos,fov,npix,ii,jj)
    function WrapAngle(a1,a2)
        x2=atan(a1,a2);
        x1=2.0*atan(0.5*sqrt(a1^2+a2^2));
        x1,x2
    end
    t,r,θ,ϕ=pos;
   
    MetricDown=gdown(a,pos);
    gtt=MetricDown[1,1];
    grr =MetricDown[2,2];
    gThetaTheta = MetricDown[3,3];
    gPhiPhi =MetricDown[4,4];
    gtPhi =MetricDown[1,4];
    grPhi=MetricDown[2,4];
    gtr=MetricDown[1,2];

    e0 =[(((-1)*grPhi^2+gPhiPhi*grr)*((-2)*grPhi*gtPhi*gtr+gPhiPhi*gtr^2+grPhi^2*gtt+grr*(gtPhi^2+(-1)*gPhiPhi*gtt))^(-1))^(1/2),(grPhi^2+(
    -1)*gPhiPhi*grr)^(-1)*((-1)*grPhi*gtPhi+gPhiPhi*gtr)*(((-1)*grPhi^2+gPhiPhi*grr)*((-2)*grPhi*gtPhi*gtr+gPhiPhi*gtr^2+grPhi^2*gtt+grr*(gtPhi^2+(-1)*gPhiPhi*gtt))^(-1))^(1/2),0,(grPhi^2+(-1)*gPhiPhi*grr)^(-1)*(grr*gtPhi+(-1)*grPhi*gtr)*(((-1)*grPhi^2+gPhiPhi*grr)
    *((-2)*grPhi*gtPhi*gtr+gPhiPhi*gtr^2+grPhi^2*gtt+grr*(gtPhi^2+(-1)*gPhiPhi*gtt))^(-1))^(1/2)]';

    e1=[0.0, -(1/sqrt(grr)), 0.0,0.0]'; 
    e2=[0.0,0.0, 1/sqrt(gThetaTheta), 0.0]'; 
    e3=[0.0,0.0,0.0,-(1/sqrt(gPhiPhi))]'; 
 
    xscr=2.0*tan(fov/2)*(ii - 0.5*(npix+1.0))/npix;
    yscr=2.0*tan(fov/2)*(jj - 0.5*(npix+1.0))/npix;
    Tx,Px=WrapAngle(xscr,yscr);
    k=1.0;
    taudot=k*(-e0+cos(Tx).*e1+sin(Tx)*cos(Px).*e2+sin(Tx)*sin(Px).*e3);
   
    momentum= MetricDown*taudot'
end
