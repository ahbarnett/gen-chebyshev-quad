% test driver for genchebquad.
% (chebfun must be in path, but only needed within genchebquad.)
% Alex Barnett 2/11/26
clear
verb=1;  % >0 gives one fig per case
tol=1e-12;

for expt=0:1        % ------------------------ loop over func family choices

  a=-1; b=1;  % std interval
  % fs = family of functions via row-vector-valued func handle
  switch expt
    case 0
      fs = @(x) x.^(0:20);  % must be able to act on on col vec of x values
    case 1
      a=0; b=1;  % interval
      pows = -0.5:0.5:20;  % non-integer powers
                           %pows = -0.55 + 0.4*(0:30);  % as in julia
      fs = @(x) x.^pows;
  end
      
  [x, w, info] = genchebquad(fs, a, b, tol, verb);   % the thing to test

  N = numel(fs((a+b)/2));  % get num funcs
  kthcol = @(y,k) y(:,k);  % get k'th col, since fs(1.0)(1) not allowed :)
  for k=1:N       % indep numerical integrals of family (quadgk needs x shape)
    I(k) = quadgk(@(x) reshape(kthcol(fs(x(:)),k),size(x)), a,b, ...
                  'abstol',1e-13,'reltol',1e-12,'maxintervalcount',1e3);
  end
  IGCQ = w'*fs(x);   % test GCQ on family
  fprintf('max |I-I_GCQ| over family = %.3g\n',max(I-IGCQ))

  if verb, figure(expt+1); set(gcf,'name',sprintf('expt %d',expt));
    subplot(2,2,1); t = linspace(a,b,1e4)';   % plot grid, must be col
    plot(t,fs(t),'-'); ysc = 2; axis([a b -ysc ysc]); title('input func set');
    subplot(2,2,2);
    plot(info.x, info.ff, '.-'); title('u o.n. funcs at dense nodes');
    subplot(2,2,3); plot([x x]', [0*w w]', 'b-'); hold on;
    plot(x,w,'b.', 'markersize',20); title('GCQ rule: w_j at each x_j');
    subplot(2,2,4);
    semilogy(1:N,abs(I-IGCQ),'.', 'markersize',10); hline(tol,'r-');
    axis([1 N 1e-16 1]); title('abs I err over func set');
  end
end
