include("../src/HausdorffBEM.jl")
using .HausdorffBEM
α = 1/3
level = 6
index = repeat([false],inner=level)
Γ = CantorDust(α,index,index)
k = 22.

Φ = Helmholtz3DscreenKernel(Γ, k)
r(x1::Float64,x2::Float64,y1::Float64,y2::Float64) = sqrt((x1-y1)^2 + (x2 - y2)^2)
Φ_nearmiss_fn(x1::Float64,x2::Float64,y1::Float64,y2::Float64)  = exp.(im*k*r(x1,x2,y1,y2))./(4π*k*r(x1,x2,y1,y2))
Φ_nearmiss = KernelBasic(Φ_nearmiss_fn)

# get different values of inner products, and compare:
Q = 6
I_ss = InnerProduct(Φ, Γ, Γ, Q)
I_nm = InnerProduct(Φ_nearmiss, Γ, Γ, Q)
I_AM =  0.002754009113161 + 0.000000004743187im
I_m = 0.002754009113161
