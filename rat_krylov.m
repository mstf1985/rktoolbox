function [V, K, H] = rat_krylov(varargin)
%RAT_KRYLOV    Rational Krylov method.
%
% Runs either the rational Krylov method for A and b with poles
% xi, or, if K and H are given instead of xi, the rational Arnoldi
% recursion imposed by (H, K) with A and b. 
% It is possible to use only real arithmetic if A and b are real
% and the poles occur in complex conjugate pairs (use CPLXSORT to
% order the poles accordingly).
%  
% Possible calls are as follows:
% 
%   [V, K, H] = rat_krylov(A, b, xi)
%   [V, K, H] = rat_krylov(A, B, b, xi)
%   [V, K, H] = rat_krylov(A, b, xi, param)
%   [V, K, H] = rat_krylov(A, B, b, xi, param)
%   [V, K, H] = rat_krylov(A, b, xi, flag)
%   [V, K, H] = rat_krylov(A, B, b, xi, flag)
%   [V, K, H] = rat_krylov(A, b, xi, flag, param)
%   [V, K, H] = rat_krylov(A, B, b, xi, flag, param)
%   [V, K, H] = rat_krylov(A, b, K, H)
%   [V, K, H] = rat_krylov(A, b, K, H, param)
%   [V, K, H] = rat_krylov(A, B, b, K, H)
%   [V, K, H] = rat_krylov(A, B, b, K, H, param)
% 
% Input arguments are:
%
%   - A, B,   N-by-N matrices;
%   - b,      N-by-1 column-vector;
%   - xi,     1-by-m row-vector, with m < N;
%   - K, H,   (m+1)-by-m upper-Hessenberg matrices;
%   - flag,   should be 'real' to trigger real arithmetic for real
%             matrices with complex-conjugate pairs of poles (all 
%             values other than 'real' will be ignored);
%   - param,  parameter structure, see below.
%
% Output arguments are:
% 
%   - V,      N-by-(m+1) matrix spanning the rational Krylov space;  
%   - K, H,   (m+1)-by-m upper-Hessenberg matrices.
%
% Using param additional features can be specified. For those that
% are not specified the default is used. Currently supported:
%
%   - param.inner_product,   function handle, default @(x, y) y'*x
%   - param.linear_solver,   function handle, default @(M, x) M\x
%  
% The rational Krylov method is described in
%
% [1] A. Ruhe. Rational Krylov: A practical algorithm for large sparse 
%     nonsymmetric matrix pencils, SIAM J. Sci. Comput., 19(5):1535--1551, 
%     1998.
%
% [2] A. Ruhe. The rational Krylov algorithm for nonsymmetric eigenvalue
%     problems. III: Complex shifts for real matrices, BIT,
%     34:165--176, 1994.


  [N, m, H, K, U, xi, mu, rho, eta, A, B, b, realopt, rerun, ...   
   inner_product, induced_norm, linear_solver] = parse_argin(varargin);

  V = zeros(N, m+1);
  
  if rerun
    rerun_krylov;
  else
    run_krylov;
  end
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Nested functions. %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
function run_krylov  
  bd = false;
  bd_tol = eps(1);
  % Starting vector.
  V(:, 1) = b/induced_norm(b);
  j = 1;
  while j <= m
    % Computing the continuation combination.
    if j == 1, U(1, 1) = 1; else
      [Q, ~] = qr(K(1:j, 1:j-1)/mu(j) - H(1:j, 1:j-1)/xi(j));      
      U(1:j, j) = Q(:, end);
    end

    % Compute new vector.
    w = V(:, 1:j)*U(1:j, j);
    w = rho(j)*(A*w) - eta(j)*(B*w);
    w = linear_solver(B/mu(j)-A/xi(j), w);
    
    % Orthogonalization.
    if isreal(xi(j)/mu(j)) || realopt == 0          
      % MGS
      for reo = 0:1
        for reo_i = 1:j
          hh(1) = inner_product(w, V(:, reo_i));
          w = w - V(:, reo_i)*hh(1);
          H(reo_i, j) = H(reo_i, j) + hh(1);
        end
      end
      H(j+1, j) = induced_norm(w);
      V(:, j+1) = w/H(j+1, j);
      
      if abs(H(j+1, j)) < bd_tol*norm(H(1:j, j)), bd = true; break; end
      
      % Setting the decomposition.
      K(1:j+1, j) = H(1:j+1, j);
      K(1:j+1, j) = rho(j)*[U(1:j, j); 0] + K(1:j+1, j)/xi(j);
      H(1:j+1, j) = eta(j)*[U(1:j, j); 0] + H(1:j+1, j)/mu(j);
    else
      V(:, j+1) = real(w);
      V(:, j+2) = imag(w);
      
      % MGS
      for j = j:j+1
        for reo = 0:1
          for reo_i = 1:j
            hh(1) = inner_product(V(:, j+1), V(:, reo_i));
            V(:, j+1)   = V(:, j+1) - V(:, reo_i)*hh(1);
            H(reo_i, j) = H(reo_i, j) + hh(1);
          end
        end
        H(j+1, j) = induced_norm(V(:, j+1));
        V(:, j+1) = V(:, j+1)/H(j+1, j);
        if abs(H(j+1, j)) < bd_tol*norm(H(1:j, j)), bd = true; break; end
      end
      
      % Setting the decomposition.
      rxi = real(1/xi(j-1)); ixi = imag(1/xi(j-1));
      cxi = [rxi ixi; -ixi rxi];

      rcnt = [real(U(1:j-1, j-1)) imag(U(1:j-1, j-1)); 0 0; 0 0];
      K(1:j+1, j-1:j) = H(1:j+1, j-1:j);
      K(1:j+1, j-1:j) = K(1:j+1, j-1:j)*cxi     + rho(j-1)*rcnt;
      H(1:j+1, j-1:j) = H(1:j+1, j-1:j)/mu(j-1) + eta(j-1)*rcnt;
    end % realopt
        
    j = j+1;             
  end % while j <= m
      
  if bd == true
    warning(['rat_krylov: ''lucky breakdown'' occured at ' ...
             'iteration ' num2str(j)]);
    V = V(:, 1:j);
    K = K(1:j, 1:j-1);
    H = H(1:j, 1:j-1);
  end
      
end % run_krylov

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function rerun_krylov
  % Starting vector.
  V(:, 1) = b;  
  j = 1;  
  while j <= m
    % Computing the continuation combination.
    if isreal(xi(j)/mu(j)) || realopt == 0
      U(1:j, j) = K(1:j,j)/mu(j) - H(1:j,j)/xi(j);
    elseif realopt == 1 || realopt == 2
      rxi = real(1/xi(j));    ixi = imag(1/xi(j));
      cxi = [rxi ixi; -ixi rxi]; iii = [1 0; 0 1];
      iii = rho(j)*iii/mu(j)-eta(j)*cxi;
      zuu = K(1:j, j:j+1)/mu(j) - H(1:j, j:j+1)*cxi;
      zuu = zuu/iii;
      U(1:j, j) = zuu(1:j, 1) + 1i*zuu(1:j, 2);      
    end
    
    % Next vector.
    w = V(:, 1:j)*U(1:j, j);
    w = rho(j)*(A*w) - eta(j)*(B*w);
    w = linear_solver(B/mu(j)-A/xi(j), w);
    
    % Orthogonalization.
    if isreal(xi(j)/mu(j)) || realopt == 0
      r = rho(j)*H(1:j+1, j) - eta(j)*K(1:j+1, j);
      % MGS simulation.
      V(:, j+1) = w - V(:, 1:j)*r(1:j, 1);
      V(:, j+1) = V(:, j+1)/r(j+1, 1);
    else
      if realopt == 1
        V(:, j+1) = real(w);
        V(:, j+2) = imag(w);
      elseif realopt == 2
        V(:, j+1) = w;

        w = V(:, 1:j)*conj(U(1:j, j));
        w = rho(j)*(A*w) - eta(j)*(B*w);
        w = linear_solver(B/mu(j)-A/conj(xi(j)), w);        
        V(:, j+2) = w;
        
        V(:, j+1:j+2) = V(:, j+1:j+2)*[0.5 -0.5i; 0.5 0.5i];
      end
      % MGS simulation.
      r = rho(j)*H(1:j+2, j:j+1) - eta(j)*K(1:j+2, j:j+1);
      r = r/iii;
      for jj = 0:1
        V(:, j+jj+1) = V(:, j+jj+1) - V(:, 1:j+jj)*r(1:j+jj, jj+1);
        V(:, j+jj+1) = V(:, j+jj+1)/r(j+jj+1, jj+1);
      end
      j = j+1;
    end    
    
    j = j+1;
  end % while j <= m
end % rerun_krylov

end % rat_krylov



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [N, m, H, K, U, xi, mu, rho, eta, ...
          A, B, b, realopt, rerun, ...
          inner_product, induced_norm, linear_solver] = parse_argin(varargin)
%PARSE_ARGIN    Process the input argument list to rat_krylov.

  msg = 'rat_krylov:parse_argin: ';
  inputtype = length(varargin{1});
  
  assert(inputtype >= 3);
  
  A = varargin{1}{1};    
  N = size(A, 1);
  assert(size(A, 2) == N, [msg 'A needs to be a square matrix']);
  
  switch inputtype
   case {3}    
    % rat_krylov(A, b, xi)    
    B  = speye(N);
    b  = varargin{1}{2};
    xi = varargin{1}{3};
    m  = size(xi, 2);
    H  = zeros(m+1, m);
    K  = zeros(m+1, m);
    inner_product = @(x, y) y'*x;
    realopt   = 0;
    rerun     = 0;      
    
   case {4}            
    if isstruct(varargin{1}{4})
      % rat_krylov(A, b, xi, param) 
      B  = speye(N);
      b  = varargin{1}{2};
      xi = varargin{1}{3};
      m  = size(xi, 2);
      H  = zeros(m+1, m);
      K  = zeros(m+1, m);
      if isfield(varargin{1}{4}, 'inner_product') 
        inner_product = varargin{1}{4}.inner_product;
      else 
        inner_product = @(x, y) y'*x; 
      end
      realopt = 0;
      rerun   = 0;
      
    elseif ischar(varargin{1}{4})
      % rat_krylov(A, b, xi, flag)
      B  = speye(N);
      b  = varargin{1}{2};
      xi = varargin{1}{3};      
      m  = size(xi, 2);
      H  = zeros(m+1, m);
      K  = zeros(m+1, m);
      inner_product = @(x, y) y'*x;
      if strcmp('real', varargin{1}{4}) realopt = 1;
      else                              realopt = 0; end
      rerun = 0;
      
    elseif size(varargin{1}{4}, 1) == 1
      % rat_krylov(A, B, b, xi)    
      B  = varargin{1}{2};
      b  = varargin{1}{3};
      xi = varargin{1}{4};
      m  = size(xi, 2);
      H  = zeros(m+1, m);
      K  = zeros(m+1, m);
      inner_product = @(x, y) y'*x;
      realopt   = 0;
      rerun     = 0;
      
    elseif size(varargin{1}{4}, 1) > 1
      % rat_krylov(A, b, K, H)    
      B  = speye(N);
      b  = varargin{1}{2};
      K  = varargin{1}{3};
      H  = varargin{1}{4};
      m  = size(H, 2);
      xi = zeros(1, m);
      realopt = 0;
      if m > 1 && any(diag(K, -2))
        realopt = 1;
      end
      inner_product = @(x, y) y'*x;
      rerun     = 1;
    else
      error([msg 'invalid number of input arguments']);      
    end
   case {5}    
    if ischar(varargin{1}{5})
      % rat_krylov(A, B, b, xi, flag)
      B  = varargin{1}{2};
      b  = varargin{1}{3};
      xi = varargin{1}{4};
      m  = size(xi, 2);
      H  = zeros(m+1, m);
      K  = zeros(m+1, m);
      if strcmp('real', varargin{1}{5}) realopt = 1;
      else                              realopt = 0; end
      inner_product = @(x, y) y'*x;
      rerun     = 0;
      
    elseif isstruct(varargin{1}{5}) && ...
          ischar(varargin{1}{4})
      % rat_krylov(A, b, xi, flag, param)
      B  = speye(N);
      b  = varargin{1}{2};
      xi = varargin{1}{3};
      m  = size(xi, 2);
      H  = zeros(m+1, m);
      K  = zeros(m+1, m);
      if isfield(varargin{1}{5}, 'inner_product') 
        inner_product = varargin{1}{5}.inner_product;
      else 
        inner_product = @(x, y) y'*x; 
      end
      if strcmp('real', varargin{1}{4}) realopt = 1;
      else                              realopt = 0; end
      rerun   = 0;      
      
    elseif isstruct(varargin{1}{5}) && ...      
          size(varargin{1}{4}, 1) == 1
      % rat_krylov(A, B, b, xi, param)  
      B  = varargin{1}{2};
      b  = varargin{1}{3};
      xi = varargin{1}{4};
      m  = size(xi, 2);
      H  = zeros(m+1, m);
      K  = zeros(m+1, m);
      if isfield(varargin{1}{5}, 'inner_product') 
        inner_product = varargin{1}{5}.inner_product;
      else 
        inner_product = @(x, y) y'*x; 
      end
      realopt = 0;
      rerun   = 0;          
      
    elseif isstruct(varargin{1}{5})
      % rat_krylov(A, b, K, H, param)      
      B  = speye(N);
      b  = varargin{1}{2};
      K  = varargin{1}{3};
      H  = varargin{1}{4};
      m  = size(H, 2);
      xi = zeros(1, m);
      realopt = 0;
      if m > 1 && any(diag(K, -2))
        realopt = 1;
      end      
      if isfield(varargin{1}{5}, 'inner_product') 
        inner_product = varargin{1}{5}.inner_product;
      else 
        inner_product = @(x, y) y'*x; 
      end
      rerun   = 1;    
      
    elseif size(varargin{1}{5}, 1) > 1
      % rat_krylov(A, B, b, K, H)
      B  = varargin{1}{2};
      b  = varargin{1}{3};
      K  = varargin{1}{4};
      H  = varargin{1}{5};
      m  = size(H, 2);
      xi = zeros(1, m);
      
      realopt = 0;
      if m > 1 && any(diag(K, -2))
        realopt = 1;
      end       
      inner_product = @(x, y) y'*x;         
      rerun     = 1;    
    else
      error([msg 'invalid number of input arguments']);
    end
   case {6}
    if isstruct(varargin{1}{6}) && ...
          ischar(varargin{1}{5})
      % rat_krylov(A, B, b, xi, flag, param)
      B  = varargin{1}{2};
      b  = varargin{1}{3};
      xi = varargin{1}{4};
      m  = size(xi, 2);
      H  = zeros(m+1, m);
      K  = zeros(m+1, m);      
      if isfield(varargin{1}{6}, 'inner_product')
        inner_product = varargin{1}{6}.inner_product;
      else 
        inner_product = @(x, y) y'*x; 
      end
      if strcmp('real', varargin{1}{5}) realopt = 1;
      else                              realopt = 0; end      
      rerun = 0;

    elseif size(varargin{1}{5}, 1) > 1
      % rat_krylov(A, B, b, K, H, param)
      B  = varargin{1}{2};
      b  = varargin{1}{3};
      K  = varargin{1}{4};
      H  = varargin{1}{5};
      m  = size(H, 2);
      xi = zeros(1, m);
      realopt = 0;
      if m > 1 && any(diag(K, -2))
        realopt = 1;
      end       
      if isfield(varargin{1}{6}, 'inner_product')
        inner_product = varargin{1}{6}.inner_product;
      else 
        inner_product = @(x, y) y'*x; 
      end
      rerun = 1;
    else
      error([msg 'invalid number of input arguments']);
    end
   otherwise 
    error([msg 'invalid number of input arguments']);
  end
  
  assert(size(B, 1) == N && size(B, 2) == N, ...
         [msg 'B needs to be a square N-by-N matrix']);  
  assert(size(b, 1) == N && size(b, 2) == 1, ...
         [msg 'b needs to be an N-by-1 vector']);  
    
  % If A, B and/or b are not real, the real option cannot be used
  % unless we are rerunning the recursion, in which case it needs
  % to be handled differently then when all A, B and b are real.
  if realopt && (~isreal(A) || ~isreal(B) || ~isreal(b))
    if rerun
      realopt = 2;
    else
      realopt = 0;
      warning(['Warning: Ignoring ''real'' option with complex' ...
               ' matrices/vectors.']);      
    end
  end

  % Check the dimensions of H, K, and xi. Check also if H is
  % upper-Hessenberg and if K is block upper-Hessenberg. The last
  % check is not fully supported.
  assert(rerun || m < N, [msg 'number of poles cannot be greater' ...
                    ' or equal to the size of the problem']);    
  assert(size(xi, 1) == 1 && size(xi, 2) == m, ...
         [msg 'xi must be a 1-by-m row vector']);          
  assert(size(H, 1) == m+1 && size(H, 2) == m, ...
         [msg 'H of inappropriate size']);
  assert(size(K, 1) == m+1 && size(K, 2) == m, ...
         [msg 'K of inappropriate size']);  
  assert(all(~any(tril(H(3:end, 1:end-1)))), ...
         [msg 'H is not upper-Hessenberg']);
  assert(all(~any(tril(K(4:end, 1:end-1)))), ...
         [msg 'K is not quasi-upper-Hessenberg']);
  
  U  = zeros(m+1, m);    
  
  % Compute the poles in the exact order they appear in (H, K).
  if rerun, xi = rat_poles(H, K); end
   
  % Cannot use real option if the poles are not ordered canonically.
  if realopt && ~canonical_cplx(xi)
    realopt = 0;
    warning(['Warning: Ignoring ''real'' option as conj(poles(j))' ...
             ' == poles(j+1) failed']);      
  end
  
  [xi, mu, rho, eta] = poles_to_moebius(xi);
  
  induced_norm  = @(x) sqrt(inner_product(x, x));
  if isstruct(varargin{1}{end}) && isfield(varargin{1}{end}, 'linear_solver')
    linear_solver = varargin{1}{end}.linear_solver;
  else
    linear_solver = @(M, x) M\x;
  end
  
end


function y = canonical_cplx(xi)
% CANONICAL_CPLX Check if the poles xi are ordered canonically.
  
  m = length(xi);
  y = 1;  
  j = 1;
  while j <= m
    if isreal(xi(j)) || isinf(xi(j))
      j = j+1;
    else
      if j == m || (j < m && xi(j+1) ~= conj(xi(j)))
        y = 0;
      end % if
      j = j+2;
    end % if
  end % while
  
end


function [xi, mu, rho, eta] = poles_to_moebius(xi)
% POLES_TO_MOEBIUS Moebius transformation with poles xi.
%
% Nonzero xi is replaced with (xi, mu) := (xi, 1) and 
% (rho, eta) := (1, 0),  and xi = 0 is replaced by
% (xi, mu) := (0, inf) and (rho, eta) := (0, 1).
  mu  = ones(size(xi));
  rho = ones(size(xi));
  eta = zeros(size(xi));
  %rho = randn(size(xi));
  %eta = randn(size(xi));
  
  mu(xi == 0) = inf;
  xi(xi == 0) = 1;
  
  eta(isinf(mu)) = 1;
  rho(isinf(mu)) = 0;  
end
