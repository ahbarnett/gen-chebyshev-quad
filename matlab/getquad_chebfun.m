function [x,w] = getquad_chebfun(f)
  % GETQUAD_CHEBFUN  spit out quadrature nodes and weights used by a chebfun
  %
  % [x,w] = getquad_chebfun(f) returns nodes x and corresponding weights w
  %  appropriate for accurate quadrature of the chebfun f.
  %
  % f may be a vector-valued chebfun (quasimatrix), in which case just the
  % first col f(:,1) is used. (The discretization is the same for all cols).
  %
  % Calling with no arguments does self-test.
  %
  % Needs Chebfun package to be in the path.

  % Mostly Dan Fortunato, expanded by Alex Barnett, 2/11/26
  if nargin==0, test_getquad_chebfun; return; end
  
  assert(numel(f) == 1)     % don't handle an array of indep chebfuns
  f = f(:,1);               % just use 1st col
  % Chebyshev grid on each piece of f, in cell array
  xc = f.points;
  if ~iscell(xc), xc = {xc}; end   % handle no-split case, ugh
  % Assemble the Chenshaw-Curtis weights and scale to the size of each piece
  wc = cell(size(xc));      % num_panels * 1
  for k = 1:numel(xc)       % loop over panels
    L = xc{k}(end) - xc{k}(1);
    wc{k} = chebtech2.quadwts(length(xc{k})).' * L/2;
  end
  % stack cell arrays to single col
  x = cat(1,xc{:,1});
  w = cat(1,wc{:,1});
end

%%%%%%
function test_getquad_chebfun
  f = chebfun(@(x) [x^x, sqrt(x)], [0 1], 'splitting', 'on');
  N = size(f,2);    % how many cols
  [x,w] = getquad_chebfun(f);
  fprintf('%d panels, %d nodes\n', numel(f(:,1).pointValues), numel(x))
  % get values at nodes (internal to a chebfun)
  ff = f.values;
  ff = reshape(cat(1, ff{:}), [], N);   % m*N array, for m fine nodes
  I = sum(ff .* w);      % vector of GCQ integrals, broadcasts w over cols
  fprintf('|our quad - chebfun quad| = %g (should match to emach)\n',...
          max(abs(I - integral(f))))
end
