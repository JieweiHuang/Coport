include("GetMetric.jl")

function ZamoProjection(a,x)
    MetricDown=gdown(a,x);
    gtt=MetricDown[1,1];
    grr =MetricDown[2,2];
    gThetaTheta = MetricDown[3,3];
    gPhiPhi =MetricDown[4,4];
    gtPhi =MetricDown[1,4];
    grPhi=MetricDown[2,4];
    gtr=MetricDown[1,2];

    e0 =[(((-1)*grPhi^2+gPhiPhi*grr)*((-2)*grPhi*gtPhi*gtr+gPhiPhi*gtr^2+grPhi^2*gtt+grr*(gtPhi^2+(-1)*gPhiPhi*gtt))^(-1))^(1/2),(grPhi^2+( 
    -1)*gPhiPhi*grr)^(-1)*((-1)*grPhi*gtPhi+gPhiPhi*gtr)*(((-1)*grPhi^2+gPhiPhi*grr)*((-2)*grPhi*gtPhi*gtr+gPhiPhi*gtr^2+grPhi^2*gtt+grr*(gtPhi^2+(-1)*gPhiPhi*gtt))^(-1))^(1/2),0,(grPhi^2+(-1)*gPhiPhi*grr)^(-1)*(grr*gtPhi+(-1)*grPhi*gtr)*(((-1)*grPhi^2+gPhiPhi*grr) 
    *((-2)*grPhi*gtPhi*gtr+gPhiPhi*gtr^2+grPhi^2*gtt+grr*(gtPhi^2+(-1)*gPhiPhi*gtt))^(-1))^(1/2)];

    e3 = [0.0,(1/sqrt(grr)), 0.0,0.0]; 
    e2 =[0.0,0.0,- 1/sqrt(gThetaTheta), 0.0]; 
    e1 = [0.0,0.0,0.0,(1/sqrt(gPhiPhi))]; 
    
  inv([e0;;e1;;e2;;e3]);
    
end


