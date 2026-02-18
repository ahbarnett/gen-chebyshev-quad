function [x, w, info] = genchebquad_dev(fs, ab, tol)
%GENCHEBQUAD_DEV   Generalized Chebyshev quadrature (GCQ) custom rule.
%
% [x,w] = genchebquad_dev(fs, ab, tol) returns nodes x and corresponding
%  weights w that approximately integrate, to relative accuracy tol, all
%  functions in the span of the function set fs, over the interval ab = [a,b].
%  fs should accept a scalar and output the set as elements of a row vector.
%  The function set can be highly linearly dependent, and weakly singular.
%
% ab may enforce breakpoints by having >2 elements (as in chebfun constructor)
%
% Example (non-integer powers):
%  [x,w] = genchebquad_dev(@(x) x.^(-0.5:0.4:20), [0 1e-12 1], 1e-10)
% still can't beat 1e-7, due to either integrals or chebfun SVD ?
%
% [x,w,info] = genchebquad_dev(fs, ab, tol) also prints diagnostics and
%  returns diagnostic info:
%   info.x, info.w = m "fine" nodes and weights implicit in chebfun approx
%   info.ff = m*N matrix of func vals of U, orthonormal basis for fs func set
%   info.inds = indices of fine nodes selected, etc
%
% Called without arguments, does self-test (using internal code, not the std test)
%
% See also: TEST_GENCHEBQUAD
%
% Notes:
%  1. This is a development implementation, and is unstable, hard to read

% Alex Barnett 2/17/26
if nargin==0, local_test_genchebquad; return; end
if nargin<3 || isempty(tol), tol=1e-10; end
verb = (nargout>2);

sl = 50;   % max num nodes per panel (typ seems to be < sl/2), and regen nodes.
if numel(ab)==2    % make quasimatrix (each col a chebfun)...
  A = chebfun(fs,ab, 'splitting','on', 'splitLength',sl);
else
  for k=1:numel(ab)-1   % forces more splits than simply numel(ab)>2 to chebfun..
    A{k} = chebfun(fs,ab([k k+1]), 'splitting','on', 'splitLength',sl);
  end
  A = join(A{:});     % splats cells to arg list, join concats chebfuns
end

% do everything by hand using more resampled nodes... (using the panel split)
[xx,ww] = getquad_chebfun(A,sl);       % regen sl nodes per panel
N = size(A,2);        % how many funcs
ff = fs(xx);          % from now we ignore the chebfun A...
ff = ff .* sqrt(ww);  % row-wise scale so l^2 norm of columns approx L^2 norms
[U,S,V] = svd(ff);
r = sum(diag(S) > S(1)*tol);   % eps-rank
if verb, fprintf('\tGCQ dev: fine nodes m=%d, eps-rank r=%d\n',numel(xx), r); end
U = U(:,1:r);          % keep o.n. cols, eps-acc (in L^2) basis for fs
[Q,R,inds] = qr(U', 'vector');   % CPQR of transpose, rank-revealing
inds = inds(1:r);    % N skeleton cols of ff' that span its range (use ID?)
inds = sort(inds);   % optional
x = xx(inds);         % convert pivot indices to real numbers x in [a,b]
info.x = xx; info.U = U; info.inds = inds; info.Q = Q; info.R = R;

%U.points{1}  % turns out a 1-pt smooth piece can become 2-pt when take SVD :(
%[x, info] = rowpivots_chebfun(U,sl);   % since no LQ factorization for chebfuns

if 1
  I = sum(U .* sqrt(ww), 1);            % use new quadr scheme not on chebfun
  Ux = U(inds,:) ./ sqrt(ww(inds));     % recover U vals at pivot nodes
  w = Ux' \ I(:);     % solve Vandermonde-transpose system for wei
  if verb, fprintf('\tw solve: cond=%.3g, rel resid nrm=%.3g, |w|_1=%.3g\n',...
                   cond(Ux), norm(Ux'*w-I)/norm(I), norm(w,1)); end


  
elseif 0         % old dead case where U was chebfun:  use o.n. sys
  Ux = U(x);             % (re!)evaluate r*r mat of o.n. col funcs at GCQ nodes
                         %Ux = info.ff(info.inds,:);  % or reuse, same
  I = integral(U)';       % row vec of true integrals of the o.n. col funcs
  %Ie = (IA*V)./diag(S)'; Ie = Ie(1:r)'  % try get I from original fs ints & SVD? unstable
  %norm(Ie-I,'inf')
  %I(1:6) = Ie(1:6);  % selective overwrite of more sing ones? didn't help
  w = Ux' \ I;     % solve Vandermonde-transpose system for wei
  %w = linsolve(Ux', I, struct('RECT',true));  % no difference
  if verb, fprintf('\tw solve: cond=%.3g, rel resid nrm=%.3g, |w|_1=%.3g\n',...
                   cond(Ux), norm(Ux'*w-I)/norm(I), norm(w,1)); end
  info.I = I; info.S = S; info.V = V;  % add to diagnostics

else    % also case where A and U were chebfuns, this failed
  a = ab(1); b = ab(end);       % get better integrals, maybe override I calc...
  N = numel(fs((a+b)/2)); kthcol = @(y,k) y(:,k);
  IA = nan(1,N);  % row
  for k=1:N
    IA(k) = quadgk(@(x) reshape(kthcol(fs(x(:)),k),size(x)), a, b,...
                   'abstol',1e-13,'reltol',1e-12,'maxintervalcount',1e4);
  end
  %IA = integral(A)'; %also good to >12 digits
  Ax = A(x);       % row-piv nodes from U above
  w = Ax' \ IA(:);     % LSQ solve rect Vandermonde-transpose system for wei
  info.IA = IA;
end



% Dev notes:
% By overriding the orthog step, setting U=A, I found info.I is not
% correct for the singular func x^-0.5, when computed as a chebfun integral
% using a standard-splitted chebfun on [0,1].
% Or something inacc in chebfun for SVD of chebfun stack w/ singular funcs?

% precond variant, no help:
%colsc = sqrt(info.w(info.inds));  % R-precond w/ fine?
%w = (Ux' * diag(colsc)) \ rhs;     % solve Vandermonde-transpose system for wei
%w = diag(colsc)*w;                % undo any R-precond

% notes on chebfun splitting criterion: I can't get splits smaller than 1e-14 :(
% The first [0 2e-14] interval only has 1 node in it! How make more?
% help chebfun.techpref
% help happinesscheck
% edit standardCheck.m <- shows how data.hscale (1/intervallen) scales tol.
% Happiness scales coeff tail by data.hscale (interval len); we don't want this.
% Tried removing hscale factor in happiness, no effect.

% Solved acc by resampling onto sl=50 nodes per chebfun panel



%%%%%%%%%%%% tester with hacks
function local_test_genchebquad
verb = 0;

tol=1e-12;  % target tol
d = 20;     % max degree

for expt = 1 %0:3   % .... loop over function families
  ab = [-1 1];   % std interval
  switch expt
    case 0       % monomials (trivial case)
      fs = @(x) x.^(0:d);        % each col elem of x produces a row
    case 1       % non-integer power set
      %ab = [0 1]; fs = @(x) x.^(-0.3:0.4:d);
      ab = [0 1e-10 1]; fs = @(x) x.^(-0.5:0.4:d);  % harder: solved by regen nodes!
    case 2       % smooth plus log-singular times smooth
      fs = @(x) [x.^(0:d), x.^(0:d).*log(abs(x-0.6))];   % stack rows
    case 3       % smooth plus nearby pole times smooth
      z0 = 0.4+1e-4i;    % pole loc
      fs = @(x) [x.^(0:d), x.^(0:d).*real(1./(x-z0)) x.^(0:d).*imag(1./(x-z0))];
  end
  fprintf('expt %d...\n', expt);
  [x, w, info] = genchebquad_dev(fs, ab, tol);   % or use _dev version
  if verb, for j=1:numel(x)
      fprintf('x_%d = %-22.17g \t w_%d = %-22.17g\n',j,x(j),j,w(j));
  end, end
  a = ab(1); b = ab(end);  % get interval
  N = numel(fs((a+b)/2));  % how many funcs?
  kthcol = @(y,k) y(:,k);  % get k'th col, since fs(1.0)(1) not allowed :(
  Is = nan(1,N);
  for k=1:N       % indep numerical integrals of family (quadgk needs x shape)
    Is(k) = quadgk(@(x) reshape(kthcol(fs(x(:)),k),size(x)), a,b, ...
                   'reltol',tol/10, 'maxintervalcount',1e4);  % this can fail too!
                                                              % reltol crucial
  end
  IGCQ = w'*fs(x);   % apply GCQ to all in family to get row vec
  fprintf('\t%d nodes: max |I_quadgk-I_GCQ| over family = %.3g\n',...
          numel(x), max(abs(Is-IGCQ)));
end              % ....
