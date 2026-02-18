function [x, w, info] = genchebquad_basic(fs, ab, tol)
%GENCHEBQUAD_BASIC   Generalized Chebyshev quadrature (GCQ) custom rule.
%
% [x,w] = genchebquad_basic(fs, ab, tol) returns nodes x and corresponding
%  weights w that approximately integrate, to relative accuracy tol, all
%  functions in the span of the function set fs, over the interval ab = [a,b].
%  fs should accept a scalar and output the set as elements of a row vector.
%  The function set can be highly linearly dependent, and weakly singular.
%
% Example (non-integer powers):
%  [x,w] = genchebquad_basic(@(x) x.^(-0.2:0.4:20), [0 1], 1e-12)
%
% Example (smooth plus log-singular times smooth):
%  fs = @(x) [x.^(0:20), x.^(0:20).*log(abs(x-0.7))];
%  [x,w] = genchebquad_basic(@(x) x.^(-0.2:0.4:20), [-1 1], 1e-12)
%
% [x,w,info] = genchebquad_basic(fs, ab, tol) also prints diagnostics and
%  returns diagnostic info:
%   info.x, info.w = m "fine" nodes and weights implicit in chebfun approx
%   info.ff = m*N matrix of func vals of U, orthonormal basis for fs func set
%   info.inds = indices of fine nodes selected, etc
%
% Called without arguments, does self-test.
%
% See also: TEST_GENCHEBQUAD
%
% Notes:
%  1. This is a basic implementation using Chebfun, making the
%     algorithm clear, but is sometimes less accurate than would be possible.)
%  2. Needs Chebfun to be in the path.

% Alex Barnett 2/17/26
if nargin==0, verb=0; test_genchebquad; return; end
if nargin<3 || isempty(tol), tol=1e-10; end
verb = (nargout>2);

A = chebfun(fs, ab, 'splitting','on', 'splitLength', 50);   % interpolate func set
% (here the degree of each chebfun panel should be grown for more accuracy...)
[U,S,V] = svd(A);                  % all N cols of U computed, where N = numel(fs)
r = sum(diag(S) > S(1)*tol);       % eps-rank
U = U(:,1:r);         % keep o.n. quasimat cols, an eps-acc (in L^2) basis for fs
[x, info] = rowpivots_chebfun(U);  % since no LQ factorization for chebfuns :(
if verb, fprintf('\tGCQ: %d fine nodes, eps-rank=%d\n',numel(info.x), r); end
Ux = U(x);            % (re)evaluate r*r mat of o.n. col funcs at GCQ nodes
I = integral(U)';     % row vec of true integrals of the o.n. funcs
w = Ux' \ I;          % solve Vandermonde-transpose for wgts
if verb, fprintf('\tw solve: cond=%.3g, rel resid nrm=%.3g, |w|_1=%.3g\n',...
                 cond(Ux), norm(Ux'*w-I)/norm(I), norm(w,1)); end
info.I = I; info.s = diag(S);      % add to diagnostics
