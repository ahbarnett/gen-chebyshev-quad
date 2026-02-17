% test driver for genchebquad.
% (chebfun must be in path, but only needed within genchebquad.)
% Chebfun 5.7.0. 
% Alex Barnett 2/11/26 - 2/16/26.
clear
verb=1;  % >0 gives one fig per case
tol=1e-12;

for expt=1        % ------------------------ loop over func family choices

  ab = [-1 1];  % std interval
  % fs = family of functions via row-vector-valued func handle
  switch expt
    case 0
      fs = @(x) x.^(0:20);  % must be able to act on on col vec of x values
    case 1
      ab = [0 1e-10 1];  % interval with helpful split to force finer chebs
                         %ab = [0 1];
                         %pows = -0.5:0.5:20;  % non-integer powers (-0.5 too bad)
      pows = -0.35 + 0.4*(0:30);  % as in julia, but quadgk fails for 1st func (power < -.5 bad)
      fs = @(x) x.^pows;
    case 2
      d = 20; x0 = 0.6;
      fs = @(x) [x.^(0:d), x.^(0:d).*log(abs(x-x0))];
  end
      
  [x, w, info] = genchebquad(fs, ab, tol, verb);   % the thing to test
  for j=1:numel(x), fprintf('x_%d = %-22.17g \t w_%d = %-22.17g\n',j,x(j),j,w(j)); end
  
  a = ab(1); b = ab(end);  % interval endpts
  N = numel(fs((a+b)/2));  % call midpt to get num funcs
  kthcol = @(y,k) y(:,k);  % get k'th col, since fs(1.0)(1) not allowed :)
  for k=1:N       % indep numerical integrals of family (quadgk needs x shape)
    I(k) = quadgk(@(x) reshape(kthcol(fs(x(:)),k),size(x)), a,b, ...
                  'reltol',tol,'maxintervalcount',1e4);
  end
  IGCQ = w'*fs(x);   % test GCQ on family
  fprintf('max |I-I_GCQ| over family = %.3g\n',max(abs(I-IGCQ)))

  if verb, figure(expt+1); clf; set(gcf,'name',sprintf('expt %d',expt));
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


% Dev notes:
% By overriding the orthog step, setting U=A, I found info.I is not
% correct for the singular func x^-0.5, when computed as a chebfun integral
% using a auto-splitted chebfun on [0,1].
% Something inacc in chebfun for SVD of chebfun stack w/ singular funcs.


