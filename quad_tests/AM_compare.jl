include("../src/HausdorffBEM.jl")
using .HausdorffBEM, Plots, SpecialFunctions
α = 0.4
level = 8
index = repeat([false],inner=level)
Γ = CantorLine(α,index)
k = 52.

Φ = Helmholtz2DscreenKernel(Γ, k)
Φ_nearmiss_fn(x::Float64, y::Float64) = im/4*besselh(0,1,k*abs(x-y))
Φ_nearmiss = KernelBasic(Φ_nearmiss_fn)

# get different values of inner products, and compare:
Q = 8
I_ss = InnerProduct(Φ, Γ, Γ, Q)
I_nm = InnerProduct(Φ_nearmiss, Γ, Γ, Q)
I_AM = 1.208425270251492e-05 + 3.814578602855435e-06im
