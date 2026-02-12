function [x, w, info] = genchebquad(fs, a, b, tol, verb)
%GENCHEBQUAD  Generalized Chebyshev quadrature (GCQ) rule construction
%
%   [x, w, info] = genchebquad(fs, a, b [, tol] [, verb])
%
%   fs   : function handle; scalar input -> row vector of function values
%   a,b  : interval endpoints
%   tol  : relative/absolute tolerance (default 1e-10)
%   verb : verbosity (default 0)
%
% Returns nodes x and weights w that integrate the function family in fs
%  on (a,b) to tolerance tol. info contains diagnostic data:
%   info.x, info.w = m nodes and weights implicit in the chebfun approximation
%   info.ff = m*N matrix of func vals of U, orthonormal basis for fs func set
%
% Needs Chebfun package to be in the path.
%
% See test_genchebquad for test/demo.

% Alex Barnett 2/11/26
if nargin < 4 || isempty(tol), tol = 1e-10; end
if nargin < 5 || isempty(verb), verb = 0; end

A = chebfun(fs,[a,b], 'splitting', 'on');  % quasimatrix (each col a chebf)
[U,S,V] = svd(A,0);    % unsure 0 does anything (all N cols of U computed)
r = sum(diag(S) > S(1)*tol);   % eps-rank
U = U(:,1:r);          % keep o.n. quasimat cols, eps-acc (in L^2) basis for fs
[x, info] = rowpivots_chebfun(U);   % since no LQ factorization for chebfuns
if verb, fprintf('GCQ:\tfine nodes m=%d, eps-rank r=%d\n',numel(info.x), r); end
Ux = U(x);             % (re!)evaluate r*r mat of o.n. col funcs at GCQ nodes
I = integral(U);       % row vec of true integrals of the o.n. col funcs
rhs = I(:);
w = Ux' \ rhs;         % solve Vandermonde transpose (Newton-Cotes) system
if verb, fprintf('\tw solve: cond=%.3g, rel resid nrm=%.3g, |w|_1=%.3g\n',cond(Ux), norm(Ux'*w-rhs)/norm(rhs), norm(w,1)); end
