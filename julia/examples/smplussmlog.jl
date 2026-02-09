include("../src/genchebquad.jl")

let
    a, b = -1.0, 1.0
    x0 = 0.6                                                            # location of log singularity
    # construct the rule...
    fs = x -> reduce(vcat, x^k .* [1, log(abs(x - x0))] for k = 0:20)   # family of 41 funcs
    x, w = genchebquad(fs, a, b, 1e-12)                             # rule for (-1,1) good to 12 digits
    # test this rule with a smooth function...
    f = x -> sin(1 + 3x)
    fp = x -> 3cos(1 + 3x)                                              # fp must be deriv of f
    println("smooth error = ", sum(w .* fp.(x)) - f(b) + f(a))
    # test this rule with a smooth times log function...
    f = x -> sin(3 * (x - x0)) * log(abs(x - x0))
    fp = x -> 3cos(3(x - x0)) * log(abs(x - x0)) + 3sinc(3(x - x0) / pi)
    println("smooth.log error = ", sum(w .* fp.(x)) - f(b) + f(a))
end

