function gdown(a,x)
t,r,θ,ϕ=x;
M=1.0;
Δ=r^2-2*M*r+a^2;
Σ=r^2+a^2*cos(θ)^2;
[[-(1-2*M*r/Σ),0,0,-2*M*r*a*sin(θ)^2/Σ];;
[0,Σ/Δ,0,0];;
[0,0,Σ,0];;
[-2*M*r*a*sin(θ)^2/Σ,0,0,(Δ+2*M*r*(r^2+a^2)/Σ)*sin(θ)^2]
]
end

function gup(a,x)
    t,r,θ,ϕ=x;
    M=1.0;
    Δ=r^2-2*M*r+a^2;
    Σ=r^2+a^2*cos(θ)^2;
    [[-((r^2+a^2)^2-a^2*sin(θ)^2*Δ)/(Σ*Δ),0,0,-2*M*r*a/(Σ*Δ)];;
    [0,Δ/Σ,0,0];;
    [0,0,1/Σ,0];;
    [-2*M*r*a/(Σ*Δ),0,0,(Δ-a^2*sin(θ)^2)/(Σ*Δ*sin(θ)^2)]
    ]
end
    