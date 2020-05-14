export InnerProduct

# -----------------  methods for quadrature ---------------#
function HausdorffMidRule(C::CantorLine,Q::Int64)
    M = getAllIndices(Q)
    x = Array{Float64}(undef,2^Q)
    w = Array{Float64}(undef,2^Q)
    for q=1:(2^Q)
        q_index = cat(C.index, reshape(M[q,:],Q),dims=1)
        C_sub = CantorLine(C.α,q_index)
        x[q] = midpoint(C_sub)
        w[q] = μ(C_sub)
    end
    (x,w)
end

function HausdorffMidRule(C::CantorDust,Q::Int64)
    M = getAllIndices(Q)
    x = Array{Float64}(undef,4^Q)
    y = Array{Float64}(undef,4^Q)
    w = Array{Float64}(undef,4^Q)
    count = 0
    for q_x=1:(2^Q)
        q_index_x = cat(C.index_x, reshape(M[q_x,:],Q), dims=1)
        for q_y=1:(2^Q)
            count = count + 1
            q_index_y = cat(C.index_y, reshape(M[q_y,:],Q), dims=1)
            C_sub = CantorDust(C.α, q_index_x, q_index_y)
            (x[count],y[count]) = midpoint(C_sub)
            w[count] = μ(C_sub)
        end
    end
    return (x,y,w)
end

function HausdorffMidRuleGalerkin(C1::CantorLine,C2::CantorLine,Q1::Int64,Q2::Int64)
    (x,w_x) = HausdorffMidRule(C1::CantorLine,Q1::Int64)
    (y,w_y) = HausdorffMidRule(C2::CantorLine,Q2::Int64)
    X = repeat(x,2^Q2)
    W_x = repeat(w_x,2^Q2)
    Y = reshape(transpose(repeat(y,1,2^Q1)),2^(Q1+Q2),1)
    W_y = reshape(transpose(repeat(w_y,1,2^Q1)),2^(Q1+Q2),1)
    return (X,Y,W_x.*W_y)
end

function HausdorffMidRuleGalerkin(C1::CantorDust,C2::CantorDust,Q1::Int64,Q2::Int64)
    (x_1, x_2, w_x) = HausdorffMidRule(C1::CantorDust,Q1::Int64)
    (y_1, y_2, w_y) = HausdorffMidRule(C2::CantorDust,Q2::Int64)
    X_1 = repeat(x_1, 4^Q2)
    X_2 = repeat(x_2, 4^Q2)
    W_x = repeat(w_x, 4^Q2)
    Y_1 = reshape(transpose(repeat(y_1,1,4^Q1)), 4^(Q1+Q2),1)
    Y_2 = reshape(transpose(repeat(y_2,1,4^Q1)), 4^(Q1+Q2),1)
    W_y = reshape(transpose(repeat(w_y,1,4^Q1)), 4^(Q1+Q2),1)
    return (X_1, X_2, Y_1, Y_2, W_x.*W_y)
end

HausdorffMidRuleGalerkin(C1::Cantor,C2::Cantor,Q::Int64) = HausdorffMidRuleGalerkin(C1::Cantor,C2::Cantor,Q::Int64,Q::Int64)

# methods for computing inner products:

function InnerProduct(K::KernelSingInfo, C1::CantorLine, C2::CantorLine, Q::Int64)
    (x,y,w) = HausdorffMidRuleGalerkin(C1,C2,Q)
    if C1==C2
        ℓ = length(C1.index)
        #println(sum(w.* K.φ.(x,y)) )
        I = sum(w.* K.φ.(x,y))  + K.J_k(ℓ)#- log(K.k/2)/(2π*2^(2ℓ)) - SingPart(C1, K.I0)/(2π)
    else
        I = sum(w.*K.Φ.(x,y))
    end
    return I
end

function InnerProduct(K::KernelSingInfo, C1::CantorDust, C2::CantorDust, Q::Int64)
    (x_1,x_2,y_1,y_2,w) = HausdorffMidRuleGalerkin(C1,C2,Q)
    if C1==C2
        ℓ = length(C1.index_x)
        #println(sum(w.* K.φ.(x,y)) )
        I = sum(w.* K.φ.(x_1,x_2,y_1,y_2))  + K.J_k(ℓ)#- log(K.k/2)/(2π*2^(2ℓ)) - SingPart(C1, K.I0)/(2π)
    else
        I = sum(w.*K.Φ.(x_1,x_2,y_1,y_2))
    end
    return I
end

function InnerProduct(K::Kernel, C1::CantorLine, C2::CantorLine, Q::Int64)
    (x,y,w) = HausdorffMidRuleGalerkin(C1, C2, Q, Q+1)
    return sum(w.*K.Φ.(x,y))
end

function InnerProduct(K::Kernel, C1::CantorDust, C2::CantorDust, Q::Int64)
    (x1,x2,y1,y2,w) = HausdorffMidRuleGalerkin(C1, C2, Q, Q+1)
    return sum(w.*K.Φ.(x1,x2,y1,y2))
end

function InnerProduct(K::Kernel, d::DiracDelta, C::CantorLine, Q::Int64)
    (x,w) = HausdorffMidRule(C,Q)
    I = sum(w.*(K.Φ).(x,d.supp))
end

function InnerProduct(C::CantorLine, f::Function, Q::Int64)
    (x,w) = HausdorffMidRule(C,Q)
    I = sum(w.*f.(x))
end

function InnerProduct(C::CantorDust, f::Function, Q::Int64)
    (x1,x2,w) = HausdorffMidRule(C,Q)
    I = sum(w.*f.(x1,x2))
end

InnerProduct(d::DiracDelta, f::Function) = f(d.supp)

InnerProduct(d::DiracDelta, f::Function, Q::Int64) = InnerProduct(d, f)
