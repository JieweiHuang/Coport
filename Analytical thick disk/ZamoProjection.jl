include("GetMetric.jl")

function ZamoProjection(a,x)
    MetricDown=gdown(a,x);
    gtt=MetricDown[1,1];
    grr =MetricDown[2,2];
    gThetaTheta = MetricDown[3,3];
    gPhiPhi =MetricDown[4,4];
    gtPhi =MetricDown[1,4];

    e0 = sqrt(-(gPhiPhi/(gtt*gPhiPhi - gtPhi*gtPhi))).*[1.0,0.0,0.0, -(gtPhi/gPhiPhi)];
    e3 = [0.0,(1/sqrt(grr)), 0.0,0.0]; 
    e2 =[0.0,0.0,-1/sqrt(gThetaTheta), 0.0]; 
    e1 = [0.0,0.0,0.0,(1/sqrt(gPhiPhi))]; 
    
  inv([e0;;e1;;e2;;e3]);
    
end


