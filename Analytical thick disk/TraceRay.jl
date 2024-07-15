include("TraceSingleRay.jl")

function TraceRay(ν0,a,pos,fov,npix)
    out=zeros(Float64,npix,npix,4);
    Threads.@threads for ii=1:npix^2
    jj=rem(ii-1,npix);
    kk=(ii-1-jj)÷npix+1;
    jj=jj+1;
    kk=npix-kk+1;
    out[jj,npix-kk+1,:].=TraceSingleRay(ν0,a,pos,fov,npix,kk,jj)
end
out
end   


    
