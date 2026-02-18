% Test driver/demo for genchebquad_basic.
% (Chebfun must be in path, but only needed within genchebquad_basic)
% Alex Barnett, 2/17/26. Tested with Chebfun 5.7.0.
% To-do: insert non-family f(x) tests with the right singularity, as Julia.

clear
verb = 1;   % 0: short text, 1: figs + more text
tol=1e-10;  % target tol
d = 20;     % max degree

for expt = 0:3   % .... loop over function families
  ab = [-1 1];   % std interval
  switch expt
    case 0       % monomials (trivial case)
      fs = @(x) x.^(0:d);        % each col elem of x produces a row
    case 1       % non-integer power set (acc poor if pow<-0.3 with basic)
      ab = [0 1];
      fs = @(x) x.^(-0.3:0.4:d);
    case 2       % smooth plus log-singular times smooth
      fs = @(x) [x.^(0:d), x.^(0:d).*log(abs(x-0.6))];   % stack rows
    case 3       % smooth plus nearby pole times smooth
      z0 = 0.4+1e-4i;    % pole loc
      fs = @(x) [x.^(0:d), x.^(0:d).*real(1./(x-z0)) x.^(0:d).*imag(1./(x-z0))];
  end
  fprintf('expt %d...\n', expt);
  [x, w, info] = genchebquad_dev(fs, ab, tol);
  if verb, for j=1:numel(x)
      fprintf('x_%d = %-22.17g \t w_%d = %-22.17g\n',j,x(j),j,w(j));
  end, end
  a = ab(1); b = ab(end);  % get interval
  N = numel(fs((a+b)/2));  % how many funcs?
  kthcol = @(y,k) y(:,k);  % get k'th col, since fs(1.0)(1) not allowed :(
  Is = nan(1,N);
  for k=1:N       % indep numerical integrals of family (quadgk needs x shape)
    Is(k) = quadgk(@(x) reshape(kthcol(fs(x(:)),k),size(x)), a,b, ...
                  'abstol',tol, 'reltol',tol/10, 'maxintervalcount',1e4);
  end
  IGCQ = w'*fs(x);   % apply GCQ to all in family to get row vec
  fprintf('\t%d nodes: max |I_quadgk-I_GCQ| over family = %.3g\n',...
          numel(x), max(abs(Is-IGCQ)));

  if verb, figure(expt+1); clf; set(gcf,'name',sprintf('expt %d',expt));
    subplot(2,2,1); t = linspace(a,b,1e4)';   % plot grid, must be col
    plot(t,fs(t),'-'); ysc = 2; axis([a b -ysc ysc]); title('input func set');
    subplot(2,2,2);
    plot(info.x, info.ff, '.-'); title('u o.n. funcs at dense nodes');
    subplot(2,2,3); plot([x x]', [0*w w]', 'b-'); hold on;
    plot(x,w,'b.', 'markersize',20); title('GCQ rule: w_j at each x_j');
    subplot(2,2,4);
    semilogy(1:N,abs(Is-IGCQ),'.', 'markersize',10); hline(tol,'r-');
    axis([1 N 1e-16 1]); title('abs I err over func set');
  end
end              % ....
