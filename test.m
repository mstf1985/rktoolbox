
n = 100; C1 = n*gallery('tridiag',n); C1(end) = C1(end)/2;
C2 = (abs(gallery('tridiag',n))+2*speye(n))/(6*n); C2(end) = C2(end)/2;
C3 = sparse(n,n); C3(n,n) = 1; F = @(z) C1 - z*C2 + C3*z/(z-1);
Nmax  = 50;
zz = linspace(4,296,50);

Sigma = 150 + 146*exp(linspace(0,2i*pi,400));  Xi = inf;
tol  = 0; R = util_nleigs(F, Sigma, Xi, tol, Nmax); 

% %%
% K = full(spdiags([AB.beta(:)./AB.xi(:),ones(AB.N,1)],-1:0,AB.N+1,AB.N));
% H = full(spdiags([AB.beta(:),AB.sigma(1:AB.N).'],-1:0,AB.N+1,AB.N));
% 
% rat = rkfunm(K,H,AB.D,0,AB);
% 
% z = pi+1i;
% v = null((z*K - H)')'; v = v/v(1);
% [ AB.b(1,z) , v(2) ]
% 
% R1 = feval(rat,z);
%R2 = AB.eval(z,50);