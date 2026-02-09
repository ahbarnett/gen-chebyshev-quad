# Generalized Chebyshev Quadrature

author: Alex H. Barnett

A scratch repo for codes for GCQ. Such codes automatically create a 1D quadrature scheme (nodes and weights) that integrate a given family of functions to a prescribed user tolerance. They are the "toy" version of generalized Gaussian quadrature (GGQ), which may be able to half the number of nodes. However, GCQ is simple, fast, and robust.

### Example

In Julia, here we build a custom quadrature on (-1,1) for smooth functions plus functions are smooth times a logarithmic singularity at ``x0=0.6`` which lies in the integration interval:

```julia-repl
    a, b = -1.0, 1.0                                                    # interval
    x0 = 0.6                                                            # location of log singularity
    # construct the rule...
    fs = x -> reduce(vcat, x^k .* [1, log(abs(x - x0))] for k = 0:20)   # family of 41 funcs
    x, w = genchebquad(fs, a, b, 1e-12)                                 # x are nodes and w weights
    f = x -> sin(1 + 3x)                                                # smooth test func
    fp = x -> 3cos(1 + 3x)                                              # fp must be deriv of f
    println("smooth error = ", sum(w .* fp.(x)) - f(b) + f(a))
    f = x -> sin(3 * (x - x0)) * log(abs(x - x0))                       # log-singular test
    fp = x -> 3cos(3(x - x0)) * log(abs(x - x0)) + 3sinc(3(x - x0) / pi)
    println("smooth.log error = ", sum(w .* fp.(x)) - f(b) + f(a))
```

This is in ``julia/examples/``. Here is a picture of the function family (top left), orthonormal set on the adaptive grid (top right), stick plot of nodes and weights (lower left), and errors on each member of the input family (lower right):

![details of Julia GCQ output, case 3 in tests](pics/gcq_julia_case3_smplussmlog.png)

See the docstring in ``julia/src`` for another example.

### Directories

 - ``julia`` for Julia code
 - ``matlab`` for MATLAB/Octave code
 - ``pics`` for output figures

### References

 - Sec. 4.3 of: J. Bremer, Z. Gimbutas, and V. Rokhlin, "A nonlinear optimization procedure for generalized Gaussian quadratures," SIAM J. Sci. Comput. 32(4), 1761--1788 (2010).
 - App. B of: D. Malhotra and A. H. Barnett, "Efficient convergent boundary integral methods for slender bodies," J. Comput. Phys. 503, 112855 (2024).

