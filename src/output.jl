# solution functions
using GeometricalPredicates

export FarField, SLP, draw

function FarField(𝚯::Array{Float64,1}, V::Basis, coeffs::Array{Complex{Float64},1}, k::Float64)
    Q = 2
    K(θ::Float64, y::Array{Float64,1}) = exp.(-im*k*cos.(θ)*y)
    FF = Array{Complex{Float64}}(undef,length(𝚯))
    for n=1:length(𝚯)
        FF[n] = 0
        for m = 1:V.DOFs
            (y,w) = HausdorffMidRule(V.elements[m],Q)
            FF[n] = FF[n] + coeffs[m]*sum(K(𝚯[n],y).*w)
        end
    end
    return FF
end

function FarField(𝚯::Array{Float64,1}, Ψ::Array{Float64,1}, V::Basis, coeffs::Array{Complex{Float64},1}, k::Float64)
    Q = 2
    K(θ::Float64, ψ::Float64, x_1::Array{Float64,1}, x_2::Array{Float64,1}) = exp.(-im*k*(cos(θ)*x_1+sin(ψ)*x_2))
    FF = Array{Complex{Float64}}(undef,length(𝚯),length(Ψ))
    @showprogress 1 "Getting far-field data: " for n=1:length(𝚯)
        for m=1:length(Ψ)
            FF[n,m] = 0
            for v = 1:V.DOFs
                (y_1,y_2,w) = HausdorffMidRule(V.elements[v],Q)
                FF[n,m] = FF[n,m] + coeffs[v]*sum(K(𝚯[n],Ψ[m],y_1,y_2).*w)
            end
        end
    end
    return FF
end

function SLP(k::Float64, V::Basis{CantorLine}, c::Array{Complex{Float64},1}, X_1::Array{Float64, 2}, X_2::Array{Float64, 2})
    ϵ = 1e-6
    Φ(x1::Float64, x2::Float64, y::Float64) = ((x1-y)^2 + x2^2) > ϵ ? im/4*besselh(0,1,k*sqrt((x1-y)^2 + x2^2)) : 0
    Q = 5
    u_scat = zero(Array{Complex{Float64}}(undef,size(X_1)))
    @showprogress 1 "Computing scattered field " for m=1:V.DOFs
        (y,w) = HausdorffMidRule(V.elements[m],Q)
        w *= c[m]
        for n=1:length(X_1)
            u_scat[n] += sum(w.*Φ.(X_1[n],X_2[n],y))
        end
    end
    return u_scat
end

function SLP(k::Float64, V::Basis{CantorDust}, c::Array{Complex{Float64},1}, X_1::Array{Float64, 2}, X_2::Array{Float64, 2}, X_3::Array{Float64, 2}, avoidPolys::Array{Polygon2D{Point2D},1}, Q=6)
    ϵ = 1e-6
    r(x1::Float64, x2::Float64, x3::Float64, y1::Float64, y2::Float64) = ((x1-y1)^2 + (x2-y2)^2 + x3^2)
    Φ(x1::Float64, x2::Float64, x3::Float64, y1::Float64, y2::Float64) = r(x1,x2,x3,y1,y2) > ϵ ?  exp(im*k*r(x1,x2,x3,y1,y2))/(4π*k*r(x1,x2,x3,y1,y2)) : 0
    u_scat = zero(Array{Complex{Float64}}(undef,size(X_1)))
    @showprogress 1 "Computing scattered field " for m=1:V.DOFs
        (y_1,y_2,w) = HausdorffMidRule(V.elements[m],Q)
        w *= c[m]
        for n=1:length(X_1)
            notInScatterer = true
            for s in avoidPolys
                if inpolygon(s,Point(X_1[n], X_2[n]))
                    notInScatterer = false
                end
            end
            if notInScatterer
                u_scat[n] += sum(w.*Φ.(X_1[n],X_2[n],X_3[n],y_1,y_2))
            end
        end
    end
    return u_scat
end

function draw(V::Basis{CantorDust}, thicken::Float64=0.05)
    ℓ = Int(round(log(4,V.DOFs)))
    α = V.elements[1].α
    S = Array{Shape}(undef,V.DOFs)
    T = Array{Polygon2D{Point2D}}(undef,V.DOFs)
    w = α^ℓ
    square(x, y) = Shape(x .+ [-w/2,w/2,w/2,-w/2], y .+ [-w/2,-w/2,w/2,w/2])
    thickSquare(x,y) = Polygon(Point(x-w/2-thicken,y-w/2-thicken),Point(x+w/2+thicken,y-w/2-thicken),Point(x+w/2+thicken,y+w/2+thicken),Point(x-w/2-thicken,y+w/2+thicken))
    for n=1:V.DOFs
        (x,y) = midpoint(V.elements[n])
        S[n] = square(x,y)
        plot!(S[n],fillcolor=:black, legend=false)
        T[n] = thickSquare(x,y)
    end
    return T
end

function draw(V::Basis{CantorLine})
    ℓ = Int(round(log(2,V.DOFs)))
    α = V.elements[1].α
    w = α^ℓ/2
    for n=1:V.DOFs
        m = midpoint(V.elements[n])
        line = [m-w,m+w]
        plot!(line,[0,0], color=:black, legend=false)
    end
end
