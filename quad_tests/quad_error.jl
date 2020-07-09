include("../src/HausdorffBEM.jl")
using .HausdorffBEM, Plots, SpecialFunctions, LaTeXStrings

function get_errors(n::Int64, α::Float64, k::Float64, Q_range::UnitRange{Int64}, Q_ref::Int64)
    if n == 2
        Γ = CantorLine(α)
        Φ = Helmholtz2DscreenKernel(Γ, k)
        Φ_nearmiss_fn(x::Float64, y::Float64) = im/4*besselh(0,1,k*abs(x-y))
        Φ_nearmiss = KernelBasic(Φ_nearmiss_fn)
    elseif n == 3
        Γ = CantorDust(α)
        Φ = Helmholtz3DscreenKernel(Γ, k)
        r(x1::Float64,x2::Float64,y1::Float64,y2::Float64) = sqrt((x1-y1)^2 + (x2 - y2)^2)
        Φ_nearmiss_fn_(x1::Float64,x2::Float64,y1::Float64,y2::Float64)  = exp.(im*k*r(x1,x2,y1,y2))./(4π*k*r(x1,x2,y1,y2))
        Φ_nearmiss = KernelBasic(Φ_nearmiss_fn_)
    end

    #get reference solution
    I = InnerProduct(Φ, Γ, Γ, Q_ref)
    #I = InnerProduct(Φ_nearmiss, Γ, Γ, Q_ref)

    q_count = 0
    errs_SS = zeros(length(Q_range))
    errs_NM = zeros(length(Q_range))
    for q = Q_range
        q_count +=1
        errs_SS[q_count] = abs(InnerProduct(Φ, Γ, Γ, q_count)-I)/abs(I)
        errs_NM[q_count] = abs(InnerProduct(Φ_nearmiss, Γ, Γ, q_count)-I)/abs(I)
    end
    #title_string = string("Quadrature convergence, k=",k,"$\\alpha=$",round(α,digits=2), " n=",n)
    # title_string = latexstring("\$k={$(Int64(k))}, n={$(n)}, \\alpha={$(round(α,digits=2))} \$")
    # plot(Q_range,[errs_SS,errs_NM], yaxis=:log,
    #     title = title_string,#"Quadrature convergence, "*L"$k=$"*string(k)*L"$\alpha=$"*string(round(α,digits=2))*L"$n=$"*string(n),
    #     label = ["Singularity subtraction " "Near miss"],
    #     legend=:outerbottom,
    #     xtickfontsize=18,ytickfontsize=18,legendfontsize=18,titlefontsize=18,
    #     lw = 3,#line width
    #     markershape = :circle
    #     )
    # xlabel!(L"Q")
    # ylabel!("Relative error")

    return errs_SS, errs_NM
end
