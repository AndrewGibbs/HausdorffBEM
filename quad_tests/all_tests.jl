include("quad_error.jl")
using MAT
Q_ref = 6
Q_range = 1:5
# k = 1.0
# α = 1/3
n = 3

# Julia's plotting capabilities are a little ropey, so I am exporting to Matlab

for k = [1.0, 5.0, 10.0]
        for α = [0.3,0.27,0.25]
                errs_SS, errs_NM = get_errors(n, α, k, Q_range, Q_ref)
                #current()
                params = Dict(
                                "k"=>k,"alpha"=>α,"n"=>n,
                                "errs_SS"=>errs_SS,
                                "errs_NM"=>errs_NM,
                                "Q"=>collect(Q_range)
                        )
                mat_name = string("n",n,
                                        ",CF",Int64(round(1/α)),
                                        ",k",Int64(k),
                                        ".mat")
                matwrite(mat_name,params)
        end
end
