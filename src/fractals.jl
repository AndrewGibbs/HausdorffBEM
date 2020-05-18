# -----------------  methods and types for Cantor sets/dust ---------------#
export CantorLine, CantorDust, Cantor, dimH

abstract type Cantor end

# following defines a Cantor set intersected with a component of
# a prefractal to that set. Following Mattila-style indexing,
# the index is an element of ${0,1}^\ell$, and refers to the
# component of the \ell'th prefractal
struct CantorLine <: Cantor
    α::Float64
    index::Array{Bool,1}
end

#if no index is given, take intersection with level 0 prefractal:
CantorLine(α::Float64) = CantorLine(α::Float64,[])

# CantorDust follows same logic as CantorLine, except now there are
# two indices:
struct CantorDust <: Cantor
    α::Float64
    index_x::Array{Bool,1}
    index_y::Array{Bool,1}
end

CantorDust(α::Float64) = CantorDust(α::Float64,[],[])

# left endpoint of a Cantor set:
L(S::CantorLine) = 0.0 + (1-S.α)*sum(S.index.*(S.α.^(0:(length(S.index)-1))))
# midpoint of Cantor set:
midpoint(S::CantorLine) = L(S) + S.α^(length(S.index))/2
#midpoint of Cantor dust:
function midpoint(Γ::CantorDust)
    (Γx,Γy) = DustDecomp(Γ)
    (midpoint(Γx),midpoint(Γy))
end
#R(S::Cantor) = L(S) + S.α^(length(S.index))
#Hausdorff measure:
μ(S::CantorLine) = 0.5^(length(S.index))
μ(S::CantorDust) = 0.5^(length(S.index_x)+length(S.index_y))

# decompose Cantor dust into the two Cantor sets whose product
# is the Cantor dust:
function DustDecomp(C::CantorDust)
    Γx = CantorLine(C.α,C.index_x)
    Γy = CantorLine(C.α,C.index_y)
    (Γx,Γy)
end

# gets the set {0,1}^Q, useful for Mattila indexing.
# indices are ordered such that the midpoints are monotonic increasing
function getAllIndices(Q::Int64)
    x = Array{Bool}(undef,2^Q,0)
    for n=1:Q
        y = Array{Bool}(undef,0)
        for m in (false,true)
            y = [y; repeat([m],2^(n-1))]
        end
        x = cat(x,repeat(y,2^(Q-n)),dims=2)
    end
    reverse(x,dims=2)
end

#get the dimension of the fractal:
dimH(C::CantorDust) = log(4)/log(1/C.α)
dimH(C::CantorLine) = log(2)/log(1/C.α)
