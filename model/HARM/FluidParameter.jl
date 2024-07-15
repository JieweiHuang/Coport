function LocalB(a,x,p,b,fvp)
    gp=gup(a,x);
    gd=gdown(a,x);
  
    u=fvp;
    k=gp*p;
    
    B0=u'*gd*b;
    k0=u'*p;

    uB=b+B0.*u;
    uk=k+k0.*u;

    b_dot_k=uB'*gd*uk;
    b_dot_b=uB'*gd*uB;
    k_dot_k=uk'*gd*uk;

    θB=acos(max(min(b_dot_k./(sqrt(abs(k_dot_k*b_dot_b))),1.0),-1.0));

    #PitchAngle
    θB=max(min(θB,pi),0.0);

    B=sqrt(abs(b_dot_b));

    B,θB
end


