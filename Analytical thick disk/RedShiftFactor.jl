include("FourVelocity.jl");
#redshift factor
function RedShiftFactor(a,x,p)
    fvp=FourVelocity(a,x)';
   -1/(fvp*p);
end