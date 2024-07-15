function FourVelocity(a,x)
    t,r,θ,ϕ=x;
    E=1.0;
    de=r^2-2*r+a^2;
    si=r^2+a^2*cos(θ)^2;
    rf=-1;
    L=sign(-1)*a*sqrt(E^2-1.0)*sin(θ)^2;
    Q=-a^2*(E^2-1)*cos(θ)^4;
    
    ut=E*(1+2*r*(r^2+a^2)/de/si)-2*a*r*L/de/si;
    ur=sign(rf)*sqrt((E^2-1)*r^4+2*r^3+2*a^2*(E^2-1)*cos(θ )^2*r^2
        +2*((a*E-L)^2+Q)*r-a^2*Q)/si;
    uphi=E*2*a*r/de/si+L*csc(θ)^2*(de-a^2*sin(θ )^2)/de/si;
    
    [ut,ur,0,uphi];
end
    