#collection of types and functions used for Hausdorff h-BEM.

module HausdorffBEM

using SpecialFunctions, ProgressMeter, Plots

export MakeBasis, getCollocationTestBasis, hBEMsolve, Basis, DiracDelta

# Dirac delta object, for use with collocation BEM
struct DiracDelta{T}
    supp::T
end

# a basis, which may consist of anything, but popular choices are
# piecewise constants or delta functions
struct Basis{T}
    elements::Array{T,1}
    DOFs::Int64
end

include("fractals.jl")
include("kernels.jl")
include("quadrature.jl")
include("output.jl")

# -----------------  methods for constructing BEM system ---------------#

#function which creates a basis of piecewise constants
function MakeBasis(C::CantorLine, ℓ::Int64)
    DOFs = 2^ℓ
    α = C.α
    h_basis = Array{CantorLine}(undef, DOFs)
    vector_indices = getAllIndices(ℓ)
    for n = 1:DOFs
        h_basis[n] = CantorLine(α, vector_indices[n, :])
    end
    return Basis{CantorLine}(h_basis, DOFs)
end

#function which creates a basis of piecewise constants
function MakeBasis(C::CantorDust, ℓ::Int64)
    DOFs = 4^ℓ
    α = C.α
    h_basis = Array{CantorDust}(undef, DOFs)
    vector_indices = getAllIndices(ℓ)
    counter = 0
    for n = 1:(2^ℓ)
        for m = 1:(2^ℓ)
            counter = counter + 1
            h_basis[counter] = CantorDust(
                α,
                vector_indices[n, :],
                vector_indices[m, :],
            )
        end
    end
    return Basis{CantorDust}(h_basis, DOFs)
end

#function which creates a basis of dirac deltas, for use in collocation
function getCollocationTestBasis(B::Basis)
    d_basis = Array{DiracDelta}(undef, B.DOFs)
    T = typeof(midpoint(B.elements[1]))
    for n = 1:B.DOFs
        d_basis[n] = DiracDelta{T}(midpoint(B.elements[n]))
    end
    return Basis(d_basis, B.DOFs)
end


# the main BEM solver. V should always be a basis of constant functions
# W can be a basis of constant functions (Galerkin) or deltas (colllocation),
# the InnerProduct function is overloaded so can interpret either.
function hBEMsolve(K::Kernel, f::Function, V::Basis, W::Basis, Q::Int64)
    #contruct the basis:
    DOFs = V.DOFs

    if typeof(V.elements[1])==CantorDust && dimH(V.elements[1])<=1
        print("\nWARNING: Hausdorff dimension of scatterer is too low to produce
        a non-zero scattered field. Setting solution to zero.\n")
        x = zeros(Complex{Float64},DOFs)
    else
        #construct the Galerkin system:
        A = Array{Complex{Float64}}(undef, DOFs, DOFs)
        b = Array{Complex{Float64}}(undef, DOFs)

        @showprogress 1 "Constructing BEM system " for n = 1:DOFs
            for m = 1:DOFs
                A[n, m] = InnerProduct(K, W.elements[n], V.elements[m], Q)
                b[n] = InnerProduct(W.elements[n], f, Q)
            end
        end

        #now get coefficients vector
        print("solving linear system...")
        x = A \ b
        print(" done.")
    end
    return x
end

#default to test and trial space the same, if only one basis is provided
hBEMsolve(K, f::Function, V::Basis, Q::Int64) = hBEMsolve(K, f, V, V, Q)

# following handles case where Function is provided instead of Kernel object
hBEMsolve(K::Function, f::Function, V::Basis, W::Basis, Q::Int64) = hBEMsolve(KernelBasic(K), f, V, W, Q)

end

#using .HausdorffBEM
