# Generalized Chebychev quadrature function. Still rather basic, not module nor HPC. 
# Barnett 1/13/24, switched to quadgk (since SGJ included seg request) 2/9/26

using FastGaussQuadrature, QuadGK, LinearAlgebra
using Printf

struct gcq_info
    segs
    xg
    wg
    U
end

"""
    x, w = genchebquad(fs, a, b[, tol]; verb=0)
    
Use generalized Chebychev quadrature (GCQ) to compute a vector of quadrature
nodes `x`, each nodes in (a,b), and a corresponding vector of weights `w`, that
integrate all functions defined by `fs` on (a,b), to relative tolerance `tol`.
`fs` should be a vector-valued function of a scalar argument, where each entry
of the vector defines a function in the set. `verb>0` gives text diagnostics.

# Example
```julia-repl
fs(x) = [x^r for r=(-1:30)/2]          # set of powers -1/2,0,1/2,1,... 
x,w = genchebquad(fs, 0,1, 1e-12)      # returns 20 nodes and weights
```

Notes:
1) If `fs` has complex outputs, the weights `w` will also be complex. They work,
   but are not very convenient. Real-valued `fs` is recommended (eg, send in Re,
   Im parts as separate functions, via splatting the tuple reim()... )
2) QuadGK.jl is applied to `fs` to choose the initial dense node set, since it
   can return its segments arrays (as of )

References:
* Sec. 4.3 of: J. Bremer, Z. Gimbutas, and V. Rokhlin, "A nonlinear optimization
    procedure for generalized Gaussian quadratures," SIAM J. Sci. Comput.
    32(4), 1761--1788 (2010).
* App. B of: D. Malhotra and A. H. Barnett, "Efficient convergent boundary integral
    methods for slender bodies," J. Comput. Phys. 503, 112855 (2024).
"""
function genchebquad(fs::Function, a, b, tol=1e-10; verb=0)
    I,E,segs,numeval = quadgk_segbuf_count(fs,Float64(a),Float64(b); atol=tol, rtol=tol)
    sort!(segs; lt = (s,t) -> s.a<t.a)        # reorder segs along real axis
    Nf = length(I)     # how many funcs
    z0,w0 = gausslegendre(14)       # gets twice GK(7,15) min order of 2*7-1=13
    # doubling of order allows integrating all products to tol
    xg = [(s.a+s.b)/2 .+ (s.b-s.a)/2*z0 for s in segs]          # dense nodes
    xg = reduce(vcat, xg)               # flatten to one vector
    wg = [(s.b-s.a)/2*w0 for s in segs]          # dense weights (real >0)
    wg = reduce(vcat, wg)               # flatten to one vector
    if verb>0; println("nfuncs=",Nf,"\ttol=",tol,"\tnsegs=",length(segs),
        "   \tm (# dense nodes) =",length(xg)); end
     # expensive fill A size m*Nf...
    A = reduce(vcat, [sqrt(wj)*fs(xj) for (xj,wj) in zip(xg,wg)]')
    S = svd(A)    # reduced
    r = sum(S.S .> tol*S.S[1])         # *** allow setting rank instead of tol
    if verb>0; println("rank (# nodes chosen) = $r"); end   #"\tU_11=",S.U[1,1])
    U = S.U[:,1:r]
    F = qr(U',ColumnNorm())            # CPQR to get nodes in 1...m
    nodeinds = F.p[1:r]                # permutation vector
    x = xg[nodeinds];                  # final subset of nodes
    # Vandermonde: vals of U-funcs at these nodes...
    V = Diagonal(1.0./sqrt.(wg[nodeinds]))*U[nodeinds,:] 
    Is = transpose(sum(Diagonal(sqrt.(wg))*U, dims=1))  # col vec of u func ints
    w = V' \ Is              # solve trans Vandermonde sys to match u integrals
    # *** fs complex not recommended; here for V complex why not transp?
    if verb>0; @printf "weights: Vandermonde^T resid nrm %.3g\n" norm(V'*w - Is); end
    info = gcq_info(segs,xg,wg,U)   # save diagnostics
    x, w[:], info      # make w col vec too. not yet sorted wrt x_j
end
