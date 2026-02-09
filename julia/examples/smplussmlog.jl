include("../src/genchebquad.jl")

let
     a, b = -1.0, 1.0                                                    # interval
    x0 = 0.6                                                            # location of log singularity
    # construct the rule...
    fs = x -> reduce(vcat, x^k .* [1, log(abs(x - x0))] for k = 0:20)   # family of 42 funcs
	println("num funcs in = ", length(fs(0.0)))
    x, w = genchebquad(fs, a, b, 1e-12)                                 # x are nodes and w weights
	println("num nodes = ", length(x))
	f = x -> sin(1 + 3x)                                                # smooth test func
    fp = x -> 3cos(1 + 3x)                                              # fp must be deriv of f
    println("smooth error = ", sum(w .* fp.(x)) - f(b) + f(a))
    f = x -> sin(3 * (x - x0)) * log(abs(x - x0))                       # log-singular test
    fp = x -> 3cos(3(x - x0)) * log(abs(x - x0)) + 3sinc(3(x - x0) / pi)
    println("smooth.log error = ", sum(w .* fp.(x)) - f(b) + f(a))
end

