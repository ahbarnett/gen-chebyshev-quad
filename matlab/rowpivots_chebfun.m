function [p, info] = rowpivots_chebfun(f)
%ROWPIVOTS_CHEBFUN  Get continuous pivots of quasimatrix (vector chebfun)
%
% p = rowpivots_chebfun(f) returns N real numbers p in [a,b], the
%  chebfun interval, that are good pivots for the stack of N=size(f,2)
%  functions in the chebfun (quasimatrix) f. This makes the N*N matrix of f
%  sampled at p as well-conditioned as possible. Pivots are found by
%  approximate row-pivoted LQ factorization of the vector chebfun f, done by
%  CPQR on the transpose of the sqrt-quadrature-weighted values matrix.
%
% [p,info] = rowpivots_chebfun(f) gives diagnostic data:
%   info.x, info.w = m nodes and weights implicit in the chebfun approximation
%   info.ff = m*N matrix of func vals
%
% Note: 1) needs Chebfun package to be in path.
%  2) this is just part of a "row-pivoted LQ" factorization of a quasimatrix, ie
%  the column-pivoted QR of the transpose of a quasimatrix, for which I
%  cannot find an implementation (it is not in Chebfun...)
  
% Alex Barnett 2/11/26

N = size(f,2);
[x,w] = getquad_chebfun(f);   % extract implied quadrature scheme (m nodes)
ff = f.values;                % all func values at all nodes, m*N
if ~iscell(ff), ff = {ff}; end       % handle no-split case, ugh
ff = reshape(cat(1, ff{:}), [], N);  % convert cell array to m*N matrix
info.ff = ff; info.x = x; info.w = w;   % diagnostic output
ff = ff .* sqrt(w);  % row-wise scale so l^2 norm of columns approx L^2 norms
[Q,R,inds] = qr(ff', 'vector');   % CPQR of transpose, rank-revealing
inds = inds(1:N);    % N skeleton cols of ff' that span its range (use ID?)
p = x(inds);         % convert pivot indices to real numbers in [a,b]
