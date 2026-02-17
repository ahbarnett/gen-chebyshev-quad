# test script for generalized Chebyshev quadratures. Barnett 1/20/24, 2/9/26
using Printf
using QuadGK
using CairoMakie
using Test

include("../src/genchebquad.jl")  # the thing to test

function errquadfuncset(x::Vector, w::Vector, fs::Function, a, b)
    # returns the vector of errors of the rule (x,w) for each function in the set.
    I, _, nev = quadgk_count(fs, a, b; rtol=1e-12, atol=1e-14)  # est true ints (rtol can fail)
    Ig = reduce(hcat, fs.(x)) * w    # test rule (x,w). reduce makes Nf*n mat (n=#nodes)
    maxerr = norm(I - Ig, Inf)
    @printf "max abs I err in fs vs quadgk (which took %d evals) = %.3g\n" nev maxerr
    I - Ig
end

verb = 1   # >0 gives one figure per case
tol = 1e-12

@testset "genchebquad" begin

    for case = 1      # ------------------------------- loop over func set choices
        @printf("\ncase %d: ", case)
        x0 = 0.0    # singularity loc, needed to access later outside conditional block
        f = fp = fpnam = nothing   # allow defined outside conditional block
        # note how fs must be anon func, since used outside conditional block...
        if case == 0
            println("basic monomials...")
            a = -1.0
            b = 1.0        # interval
            fs = x -> [x^k for k = 0:20]       # each family member is output vec cmpnt
            k = 5; f = x -> sin(1 + k*x)    # k sets osc freq
            fp = x -> k*cos(1 + k*x)   # test func fp, must be deriv of f
            fpnam = "entire"                       # name of test func
        elseif case == 1
            println("non-integer power set...")
            a = 0.0
            b = 1.0
            # set of not-very-rational integrable powers >-1 ...
            #worstpow = -0.55; powstep = 0.4
            worstpow = -0.5; powstep = 0.5 
            fs = x -> [x^r for r = worstpow .+ powstep * (0:30)]
            # note: fails for min r<-0.6. not sure why - try arb prec?
            # *** to-do: devise and test integrand from this power set
        elseif case == 2
            println("poly plus nearby complex sqrt times poly...")
            a = -1.0
            b = 1.0
            p = 20                  # max degree
            z0 = 0.3 + 1e-3im    # sing loc near (a,b)
            #z0 = -1.001          # or, keeping it real, and sqrt away from its cut :)
            #fs = x-> [[x^k for k=0:p]; [real(x^k/sqrt(x-z0)) for k=0:p];
            #    [imag(x^k/sqrt(x-z0)) for k=0:p]]        # separate Re, Im parts
            # neater way but ordering interleaved...    
            #fs = x-> reduce(vcat, x^k.*[1, real(1/sqrt(x-z0)), imag(1/sqrt(x-z0))] for k=0:p)
            # even neater way using splat from tuple to Vector...
            fs = x -> reduce(vcat, x^k .* [1, reim(1 / sqrt(x - z0))...] for k = 0:p)
            #f(x) = sqrt(x-z0); fp(x) = 0.5/sqrt(x-z0)  # too easy
            f = x -> sin(1 + 3x) * sqrt(x - z0)
            fp = x -> 0.5 / sqrt(x - z0) * (sin(1 + 3x) + 6(x - z0) * cos(1 + 3x))
            fpnam = "anal times 1/sqrt"
        elseif case == 3
            println("poly plus log|x-x0| times poly (x0 in interval)...")
            a = -1.0
            b = 1.0
            p = 20   # max degree
            x0 = 0.6            # sing loc in (a,b)
            fs = x -> reduce(vcat, x^k .* [1, log(abs(x - x0))] for k = 0:p)
            #f(x) = (x-x0)*log(abs(x-x0)); fp(x) = 1 + log(abs(x-x0))   # too easy
            # note sinc is sin(pi.x)/(pi.x) in Julia...
            f = x -> sin(3 * (x - x0)) * log(abs(x - x0))
            fp = x -> 3cos(3(x - x0)) * log(abs(x - x0)) + 3sinc(3(x - x0) / pi)
            fpnam = "anal + log|x-x0|.anal"
        elseif case == 4
            println("poly/(x-x0), for real x0 near interval...")
            a = -1.0
            b = 1.0
            p = 16
            x0 = 1.0 + 1e-3
            fs = x -> reduce(vcat, x^k / (x - x0) for k = 0:p)
            f = x -> log(abs(x - x0))
            fp = x -> 1 / (x - x0)
            fpnam = "anal/(x-x0)"
       elseif case == 5
            println("non-integer sing power set not at origin...")
            a = 1.0       # x[1]=1+emach but fails to get acc for worstpow=-.5
            b = 2.0
            worstpow = -0.5; powstep = 0.5 
            fs = x -> [(x-a)^r for r = worstpow .+ powstep * (0:30)]
        end
        # *** to add continuous ranges of powers or singularity locations (Dhairya)

        @time global x, w, i = genchebquad(fs, a, b, tol; verb=1)   # build a rule
        if verb>0    # show nodes, wei as 2 cols
            [@printf("x_%-3d= %-22.17g w_%-3d= %-22.17g\n",j,x[j],j,w[j]) for j in 1:length(x)]
        end
        # display(fs.(x)) # any infs?
        Ierrs = errquadfuncset(x, w, fs, a, b)             # check the rule (on family)
        pass = norm(Ierrs, Inf) < 10tol                  # unif err over funcs in family
        pass || println(Ierrs)
        @test pass
        if !isnothing(f)                               # check rule on new integrand
            Ierr = abs(sum(w .* fp.(x)) - f(b) + f(a))
            Irelerr = Ierr / abs(f(b) - f(a))
            @printf "Rel err of rule on %s integrand: %.3g\n" fpnam Irelerr
            @test Irelerr < 10tol
        end

        if verb>0              # plot stuff from gcq_info struct i
            fig = Figure()
            ax1 = Axis(fig[1, 1], title=@sprintf("Case %d: input func set", case))
            Nf = length(fs(a))   # num input funcs
            t = range(a, b, 1000)
            F = reduce(hcat, fs.(t))  # eval all fs, to plot
            for fj in eachrow(F)
                lines!(t, real.(fj))
            end
            ax2 = Axis(fig[1, 2], title="u o.n. funcs at dense nodes")
            m, r = size(i.U)
            for j = 1:r             # cancel out sqrt-w factors to view raw u funcs
                scatterlines!(i.xg, real.(i.U[:, j]) ./ sqrt.(i.wg), markersize=5)
            end
            ax3 = Axis(fig[2, 1], title="rule: w_j at each x_j")
            scatter!(x, w)
            linesegments!(kron(x, [1; 1]), kron(w, [0; 1])) # stick plot
            ax4 = Axis(fig[2, 2], yscale=log10, limits=(nothing, (1e-16, 1)),
                title="abs I err over func set")
            scatter!(abs.(Ierrs))
            lines!([0, Nf], [tol, tol], color=:red, label="tol")
            axislegend()
            display(fig)
        end

    end                        # -------------------------------------------
end   # testset
