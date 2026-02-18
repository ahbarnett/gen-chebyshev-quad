function [x,w] = getquad_chebfun(f, n)
% GETQUAD_CHEBFUN  spit out quadrature nodes and weights implied by a chebfun
%
% [x,w] = getquad_chebfun(f) returns nodes x and corresponding weights w
%  appropriate for accurate quadrature of the chebfun f (assumes type-2 nodes).
%
% f may be a vector-valued chebfun (quasimatrix), in which case just the
% first col f(:,1) is used. (Chebfun discretization is the same for all cols.)
%
% [x,w] = getquad_chebfun(f,n) regenerates n new type-1 nodes per panel,
%  again returning nodes x and weights w, using only chebfun splitting from f.
%  n should be <1e3 or so.
%
% Calling with no arguments does self-test.
%
% Needs Chebfun package to be in the path.

% Original version by Dan Fortunato, expanded by Alex Barnett, 2/11/26
  if nargin==0, test_getquad_chebfun; return; end
  if nargin<2, n=[]; end
  newnodes = ~isempty(n);
  if newnodes, [gx,gw] = gauss(n); end   % precompute G-L
  
  assert(numel(f) == 1)     % don't handle an array of indep chebfuns
  f = f(:,1);               % just use 1st col
  % Chebyshev grid on each piece of f, in cell array
  xc = f.points;            % cell array of type-2 nodes
  if ~iscell(xc), xc = {xc}; end   % handle no-split case, ugh
  % Assemble the weights and scale to the size of each piece
  wc = cell(size(xc));      % num_panels * 1
  for k = 1:numel(xc)       % loop over panels
    L = f.domain(k+1)-f.domain(k); mid = (f.domain(k+1)+f.domain(k))/2;
    if newnodes
      xc{k} = mid + gx*L/2;   % overwrite node cells
      wc{k} = gw(:) * L/2;    % their weights
    else
      wc{k} = chebtech2.quadwts(length(xc{k})).' * L/2;  % CC weights
    end
  end
  % stack cell arrays to single col
  x = cat(1,xc{:,1});
  w = cat(1,wc{:,1});
end

%%%%%%
function test_getquad_chebfun
  fs = @(x) [x.^x, sqrt(x)];
  f = chebfun(fs, [0 1], 'splitting', 'on', 'splitLength', 50);
  N = size(f,2);    % how many cols
  np = numel(f(:,1).pointValues);
  
  [x,w] = getquad_chebfun(f);
  fprintf('orig chebfun: %d panels, m=%d nodes\n', np, numel(x))
  % get values at nodes (internal to a chebfun)
  ff = f.values;
  ff = reshape(cat(1, ff{:}), [], N);   % m*N array, for m fine nodes
  I = sum(ff .* w);      % vector of GCQ integrals, broadcasts w over cols
  fprintf('|our quad - chebfun quad| = %g (should match to emach)\n',...
          max(abs(I - integral(f))))
  
  n = 50;  % add regen-nodes test...
  [x,w] = getquad_chebfun(f,n);
  fprintf('regen mode (n=%d): now m=%d nodes\n', n, numel(x))
  ff = fs(x);      % expensive, eval at nodes, m*N
  I = sum(ff .* w);      % vector of GCQ integrals, broadcasts w over cols
  fprintf('|our quad - chebfun quad| = %g (should match to high acc)\n',...
          max(abs(I - integral(f))))
end
