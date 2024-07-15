function GetRadiationParameter(ne,Te,B,θB,ν0,g)
    #Physical constant
    ee=4.80320425e-10;
    h=6.62606885e-27;
    k=1.3806488e-16;
    c=2.99792458e10;
    me=9.1093829e-28;
    MBH=1.989e33*6.5e9;
    G=6.674e-8;
    ν=ν0/g;
    θe=k*Te/(me*c^2);
   
    #Fitting function
    function I_I(x)
        I_I= 2.5651 * (1.0 + 1.92 * x^(-1/3) + 0.9977 * x^(-2/3)) * exp(-1.8899 * x^(1/3));
    end
    
    function I_Q(x)
        I_Q= 2.5651 * (1.0 + 0.932 * x^(-1/3) + 0.4998 * x^(-2/3)) * exp(-1.8899 * x^(1/3));
    end
    
    function I_V(x)
        I_V= (1.8138 / x + 3.423 * x^(-2/3) + 0.02955 * x^(-0.5) + 2.0377 * x^(-1/3)) * exp(-1.8899 * x^(1/3));
    end

    #Bessel function approximations:
    function K_0(x)
        -log(0.5*x)-0.5772; 
    end

    function K_1(x)
        1.0/x; 
    end 

    function K_2(x)
        2.0/(x^2); 
    end



    function j_I(theta_e, n_e, nu, B, theta_B)
        nu_c = 3.0 * ee * B * sin(theta_B)* theta_e* theta_e / (4.0 * pi * me * c) ;
        x = nu / nu_c;
        j_I= n_e * ee * ee * nu / 2.0 / sqrt(3.0) / c / theta_e / theta_e *I_I(x);
    end
    
    
        function j_Q(theta_e, n_e, nu, B, theta_B)
         nu_c = 3.0* ee * B * sin(theta_B) / (4.0 * pi * me * c) * theta_e * theta_e;
         x = nu / nu_c;
        j_Q= n_e * ee * ee * nu / (2.0 * sqrt(3.0) * c * theta_e * theta_e) * I_Q(x);
        end
    
        function j_V(theta_e, n_e, nu, B, theta_B)
         nu_c = 3.0*ee * B * sin(theta_B) / (4.0 * pi * me * c) * theta_e * theta_e;
         x = nu / nu_c;
        j_V= 2.0* n_e * ee * ee * nu / tan(theta_B) / (3.0 * sqrt(3.0) * c * theta_e * theta_e * theta_e) * I_V(x);
        end
    
    function PlanckFunction(ν,T)
        B_nu=2.0*h*ν^3/(c^2)/(exp(h*ν/(k*T))-1.0);
    end
    function absorbtion(j_nu,nu,Te)
        B_nu=PlanckFunction(nu,Te); 
        alpnu=j_nu / B_nu;
    end
    
    
    
    
        function DeltaJ_5(X)
        D5= 0.4379 * log(1.0+0.001858*X^(1.503));
        end
    
        function f_m(X)
        fm= 2.011 * exp(-X^1.035 / 4.7) - cos(X * 0.5) * exp(-X^1.2 / 2.73) - 0.011 * exp(-X / 47.2) + (0.011 * exp(-X/47.2) - 2^(-1/3)*3^(-23/6) * 10000.0* pi * X^(-8/3)) * 0.5 * (1.0 + tanh(10.0*log(X / 120.0)));
        end
    
    
        function rho_Q(theta_e, n_e, nu, B, theta_B)
            
         c1      = 4.0 * pi * n_e * ee * ee / me;
            
         omega   = ee * B / me / c;
          X       = theta_e * sqrt(sqrt(2.0) * sin(theta_B) * (1000.0*omega/ 2.0/ pi / nu));
         Thetaer = 1.0 / theta_e;
    
        rQ= 2.0 * pi * nu / 2.0 / c * c1 * omega * omega / (2.0*pi*nu)^4 * f_m(X) * (K_1(Thetaer) / K_2(Thetaer) + 6.0* theta_e) * sin(theta_B) * sin(theta_B);
    
        end
    
    
        function  rho_V(theta_e, n_e, nu, B, theta_B)
         c1      = 4.0 * pi * n_e * ee * ee / me;
         omega   = ee * B / me / c;
         X       = theta_e * sqrt(sqrt(2.0) * sin(theta_B) * (1000.0 * omega / 2.0 / pi / nu));
         Thetaer = 1.0 / theta_e;
    
        rV= 2.0 * pi * nu / c * c1 * omega / (2.0 * pi * nu)^3 * (K_0(Thetaer) - DeltaJ_5(X)) / K_2(Thetaer) * cos(theta_B);
        end

jI = j_I(θe,ne,ν, B,θB);
jQ = j_Q(θe,ne,ν, B, θB);
jU = 0.0;
jV = j_V(θe,ne,ν, B, θB);

rU=0.0;

rQ =rho_Q(θe,ne,ν,B,θB);
rV =rho_V(θe,ne,ν,B,θB);


aI =absorbtion(jI,ν,Te);
aQ =absorbtion(jQ,ν,Te);
aU =0.0;
aV =absorbtion(jV,ν,Te);

#Affine parameter scale

rg=G*MBH/(c^2);

rg.*[(g^2).*[jI,jQ,jU,jV];[aI,aQ,aU,aV,rQ,rU,rV]./g];

end
