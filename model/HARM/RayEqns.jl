function RayEqns(a,x,p)
    pt,pr,pθ,pϕ=p;
    t,r,θ,ϕ=x;
    [(1/2)*(4*pr*r*(r^2+a^2*cos(θ)^2)^(-1)+(-2)*pt*( 
  r^2+a^2*cos(θ)^2)^(-1)*(r*(2+r)+a^2*cos(θ)^2)),( 
  1/2)*(2*a*pϕ*(r^2+a^2*cos(θ)^2)^(-1)+4*pt*r*(r^2+ 
  a^2*cos(θ)^2)^(-1)+2*pr*(a^2+((-2)+r)*r)*(r^2+a^2* 
  cos(θ)^2)^(-1)),pθ*(r^2+a^2*cos(θ)^2)^(-1),( 
  1/2)*(2*a*pr*(r^2+a^2*cos(θ)^2)^(-1)+2*pϕ*(r^2+ 
  a^2*cos(θ)^2)^(-1)*csc(θ)^2),(1/2)*(4*a*pϕ*pr* 
  r*(r^2+a^2*cos(θ)^2)^(-2)+2*pθ^2*r*(r^2+a^2* 
  cos(θ)^2)^(-2)+8*pr*pt*r^2*(r^2+a^2*cos(θ)^2) 
  ^(-2)+2*pr^2*r*(a^2+((-2)+r)*r)*(r^2+a^2*cos(θ)^2) 
  ^(-2)+(-4)*pr*pt*(r^2+a^2*cos(θ)^2)^(-1)+(-1)* 
  pr^2*((-2)+2*r)*(r^2+a^2*cos(θ)^2)^(-1)+pt^2*(2+2* 
  r)*(r^2+a^2*cos(θ)^2)^(-1)+(-2)*pt^2*r*(r^2+a^2* 
  cos(θ)^2)^(-2)*(r*(2+r)+a^2*cos(θ)^2)+2*pϕ^2*r* 
  (r^2+a^2*cos(θ)^2)^(-2)*csc(θ)^2),(1/2)*((-2)* 
  a^2*pϕ^2*(r^2+a^2*cos(θ)^2)^(-2)*cot(θ)+2* 
  pϕ^2*(r^2+a^2*cos(θ)^2)^(-1)*cot(θ)*csc(θ) 
  ^2+(-4)*a^3*pϕ*pr*cos(θ)*(r^2+a^2*cos(θ)^2)^( 
  -2)*sin(θ)+(-2)*a^2*pθ^2*cos(θ)*(r^2+a^2* 
  cos(θ)^2)^(-2)*sin(θ)+(-8)*a^2*pr*pt*r*cos(θ) 
  *(r^2+a^2*cos(θ)^2)^(-2)*sin(θ)+(-2)*a^2*pr^2* 
  (a^2+((-2)+r)*r)*cos(θ)*(r^2+a^2*cos(θ)^2)^(-2)* 
  sin(θ)+(-2)*a^2*pt^2*cos(θ)*(r^2+a^2*cos(θ) 
  ^2)^(-1)*sin(θ)+2*a^2*pt^2*cos(θ)*(r^2+a^2* 
  cos(θ)^2)^(-2)*(r*(2+r)+a^2*cos(θ)^2)*sin(θ)) 
  ];
end

