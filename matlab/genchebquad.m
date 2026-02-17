function [x, w, info] = genchebquad(fs, ab, tol, verb)
%GENCHEBQUAD  Generalized Chebyshev quadrature (GCQ) rule construction
%
%   [x, w, info] = genchebquad(fs, ab, [, tol] [, verb])
%
%   fs   : function handle; scalar input -> row vector of function values
%   ab : interval endpoints [a,b], or if length>2, forces separate
%           chebfuns which are joined (useful for strong singularities).
%   tol  : relative/absolute tolerance (default 1e-10)
%   verb : verbosity (default 0)
%
% Returns nodes x and weights w that integrate the function family in fs
%  on (a,b) to tolerance tol. info contains diagnostic data including:
%   info.x, info.w = m "fine" nodes and weights implicit in chebfun approx
%   info.ff = m*N matrix of func vals of U, orthonormal basis for fs func set
%   info.inds = indices of fine nodes selected, etc
%
% Needs Chebfun package to be in the path.
%
% Simple example (non-integer powers):
%  [x,w,~] = genchebquad(@(x) x.^(0.3:0.4:20), [0 1], 1e-12, 1)
%
% Example of endpoint singularity needing extra help:
%  [x,w,~] = genchebquad(@(x) x.^(-0.5:0.5:20), [0 1e-10 1], 1e-12, 1)
%
% See test_genchebquad for tests/demo.

% Alex Barnett 2/15/26
if nargin < 3 || isempty(tol), tol = 1e-10; end
if nargin < 4 || isempty(verb), verb = 0; end

sl = 160; %50;   % max num nodes per panel (typ seems to be about half this).
if numel(ab)==2    % make quasimatrix (each col a chebfun)..,
  A = chebfun(fs,ab, 'splitting','on', 'splitLength',sl);
else
  for k=1:numel(ab)-1
    A{k} = chebfun(fs,ab([k k+1]), 'splitting','on', 'splitLength',sl);
  end
  A = join(A{:});     % splats cells to arg list, concats chebfuns
end

[U,S,V] = svd(A);    % unsure 0 does anything (all N cols of U computed)
r = sum(diag(S) > S(1)*tol);   % eps-rank
U = U(:,1:r);          % keep o.n. quasimat cols, eps-acc (in L^2) basis for fs
%U.points{1}  % turns out a 1-pt smooth piece can become 2-pt when take SVD :(

[x, info] = rowpivots_chebfun(U);   % since no LQ factorization for chebfuns
                                    %x(x==min(x)) = 1e-27;  % debug, wasn't it
                                    %x = sort(x);           % ascending

if verb, fprintf('GCQ:\tfine nodes m=%d, eps-rank r=%d\n',numel(info.x), r); end

if 1         % use o.n. sys
  Ux = U(x);             % (re!)evaluate r*r mat of o.n. col funcs at GCQ nodes
                         %Ux = info.ff(info.inds,:);  % or reuse, same
  I = integral(U)';       % row vec of true integrals of the o.n. col funcs
  %Ie = (IA*V)./diag(S)'; Ie = Ie(1:r)'  % get I from original fs ints & SVD
  %norm(Ie-I,'inf')
  %I(1:6) = Ie(1:6);  % selective overwrite
  w = Ux' \ I;     % solve Vandermonde-transpose system for wei
                    %w = Ux' \ I;     % solve Vandermonde-transpose system for wei
                    %w = linsolve(Ux', I, struct('RECT',true));
  if verb, fprintf('\tw solve: cond=%.3g, rel resid nrm=%.3g, |w|_1=%.3g\n',cond(Ux), norm(Ux'*w-I)/norm(I), norm(w,1)); end
  info.I = I; info.S = S; info.V = V;  % add to diagnostics
else    % failed
  a = ab(1); b = ab(end);       % better integrals, maybe override I calc...
  N = numel(fs((a+b)/2)); kthcol = @(y,k) y(:,k);
  IA = nan(1,N);  % row
  for k=1:N
    IA(k) = quadgk(@(x) reshape(kthcol(fs(x(:)),k),size(x)), a, b,...
                   'abstol',1e-13,'reltol',1e-12,'maxintervalcount',1e4);
  end
  %IA = integral(A)'; %also good to >12 digits
  Ax = A(x);
  w = Ax' \ IA(:);     % solve rect Vandermonde-transpose system for wei
  info.IA = IA;
end




% precond variant, no help:
%colsc = sqrt(info.w(info.inds));  % R-precond w/ fine?
%w = (Ux' * diag(colsc)) \ rhs;     % solve Vandermonde-transpose system for wei
%w = diag(colsc)*w;                % undo any R-precond

% notes on chebfun splitting criterion: I can't get splits smaller than 1e-14 :(
% help chebfun.techpref
% help happinesscheck
% edit standardCheck.m <- shows how data.hscale (1/intervallen) scales tol.
% Happiness scales coeff tail by data.hscale (interval len); we don't want this.
% to stack two chebfuns as one: join.


