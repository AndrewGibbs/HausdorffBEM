abstract type Kernel end

struct KernelBasic <: Kernel
    Φ::Function
end

export KernelBasic, KernelSingInfo, Kernel, Helmholtz2DscreenKernel
export Helmholtz3DscreenKernel

struct KernelSingInfo  <: Kernel
    Φ::Function #the default value for the Kernel
    φ::Function #the kernel minus the singular part
    J_k::Function
end

function Helmholtz2DscreenKernel(C::CantorLine, k::Float64)
    Φ(x::Float64, y::Float64) = im/4*besselh(0,1,k*abs(x-y))
    small_threshold = 10^-14
    φ(x::Float64, y::Float64) = abs(x-y)>small_threshold ? Φ(x,y) + log(k*abs(x-y)/2)/(2π) : im/4 - Base.MathConstants.eulergamma/(2*π)

    # now create I_0, which is the integral of log(|x-y|) over Γ×Γ
    α = C.α
    Γ1 = CantorLine(α, [0])
    Γ2 = CantorLine(α, [1])
    Qbig = 10 # 14 is the highest my laptop can handle without laptop crashing
    (x,y,w) = HausdorffMidRuleGalerkin(Γ1,Γ2,Qbig,Qbig)
    I_0 = (1/0.5)*(0.5*log(α) +2*sum([w[i]*log(abs(x[i]-y[i])) for i in 1:length(x)]))

    #now create J_k as a function, which is the integral of Φ-φ over Γ_ℓ^j×Γ_ℓ^j

    J_k(ℓ::Int64) =  - log(k/2)/(2π*2^(2ℓ)) - (2.0)^(-2*ℓ)*(I_0 + ℓ*log(α))/(2π)
    return KernelSingInfo(Φ,φ,J_k)
end

function Helmholtz3DscreenKernel(C::CantorDust, k::Float64)
    r(x1::Float64,x2::Float64,y1::Float64,y2::Float64) = sqrt((x1-y1)^2 + (x2 - y2)^2)
    Φ(x1::Float64,x2::Float64,y1::Float64,y2::Float64)  = exp.(im*k*r(x1,x2,y1,y2))./(4π*k*r(x1,x2,y1,y2))
    small_threshold = 10^-14
    φ(x1::Float64,x2::Float64,y1::Float64,y2::Float64) = r(x1,x2,y1,y2)>small_threshold ? Φ(x1,x2,y1,y2) - 1 ./(4π*k*r(x1,x2,y1,y2)) : im/(4π)
    Qbig = 6
    α = C.α
    I_0 = 0.0
    Γ_nextLevel = Array{CantorDust}(undef,4)
    Γ_nextLevel[1] = CantorDust(α,[0],[0])
    Γ_nextLevel[2] = CantorDust(α,[1],[0])
    Γ_nextLevel[3] = CantorDust(α,[0],[1])
    Γ_nextLevel[4] = CantorDust(α,[1],[1])
    for n = 1:4
        for m = (n+1):4
            #if n != m
            (x_1,x_2,y_1,y_2,w) = HausdorffMidRuleGalerkin(Γ_nextLevel[n], Γ_nextLevel[m], Qbig, Qbig)
            I_0 += sum(w./(4π*r.(x_1,x_2,y_1,y_2)) )
            #end
        end
    end
    I_0 *= 2*4α/(4α-1)

    J_k(ℓ::Int64) = I_0/(k*4^(2ℓ)*α^ℓ)
    return KernelSingInfo(Φ,φ,J_k)
end
