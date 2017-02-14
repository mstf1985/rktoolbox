function [V, K, H, param] = rat_krylov(varargin)
% RAT_KRYLOV    Rational Krylov method.
%
% This function is at the core of the Rational Krylov Toolbox,
% released as part of [1]. There are three main
% functionalities provided by rat_krylov:
%   (I)   run the rational Arnoldi algorithm [4, 5];
%   (II)  extend a rational Arnoldi decomposition (RAD);
%   (III) rerun recursion given by an RAD [2].
%
% Each of the functionalities can be invoked for:
%   (i)   a single matrix A;
%   (ii)  a matrix pencil (A, B);
%   (iii) a matrix pencil (A, B) represented by a structure AB.
% Lastly, a structure param may be provided with additional
% options detailed below. Instead of param the string 'real' may be
% provided, which is equivalent to providing the structure param
% with only the field 'real' being set to 1.
%
% We now list the possible calls. The output is always of the form
% [V, K, H, param], and is therefore excluded from the following list.
%
%   rat_krylov(A, b, xi)                 % (I) - (i)
%   rat_krylov(A, b, xi, param)          % (I) - (i)
%   rat_krylov(A, V, K, H, xi)           % (II) - (i)
%   rat_krylov(A, V, K, H, xi, param)    % (II) - (i)
%   rat_krylov(A, V, K, H)               % (III) - (i)
%   rat_krylov(A, V, K, H, param)        % (III) - (i)
%
%   rat_krylov(A, B, b, xi)              % (I) - (ii)
%   rat_krylov(A, B, b, xi, param)       % (I) - (ii)
%   rat_krylov(A, B, V, K, H, xi)        % (II) - (ii)
%   rat_krylov(A, B, V, K, H, xi, param) % (II) - (ii)
%   rat_krylov(A, B, V, K, H)            % (III) - (ii)
%   rat_krylov(A, B, V, K, H, param)     % (III) - (ii)
%
%   rat_krylov(AB, b, xi)                % (I) - (iii)
%   rat_krylov(AB, b, xi, param)         % (I) - (iii)
%   rat_krylov(AB, V, K, H, xi)          % (II) - (iii)
%   rat_krylov(AB, V, K, H, xi, param)   % (II) - (iii)
%   rat_krylov(AB, V, K, H)              % (III) - (iii)
%   rat_krylov(AB, V, K, H, param)       % (III) - (iii)
%
% The input arguments are:
%
%   - A, B,   N-by-N matrices;
%   - AB,     struct representing an N-by-N pencil, see below;
%   - b,      N-by-n, here n = 1, column-vector;
%   - V,      N-by-n matrix;
%   - xi,     1-by-m row-vector, with m + n < N;
%   - K, H,   upper-Hessenberg matrices with one row more then columns;
%   - param,  struct, see below.
%
% The fields of the param structure are listed below. The names are
% self-explanatory, and if not provided, the default is being used.
%
% param.orth,          'MGS' (default), 'CGS';
% param.reorth,        0, 1 (default);
% param.inner_product, function handle of the form
%                      @(x, y) y'*x; (default)
% param.real,          0 (default), 1, 2;
%                      whether to use the real version [3] of
%                      rational Arnoldi or not;
%                      option 2 is used only for (III) with
%                      complex data on a quasi-RAD;
% param.refinement,    0 (default), 1;
%                      whether to do one step of iterative
%                      refinement for the linear system solves or
%                      not;
% param.waitbar,       0 (default), 1;
%                      whether to include a waitbar showing the
%                      progress of rat_krylov or not;
% param.continuation,        'ruhe' (default), 'almost-last', 'last',
%                            'near-optimal'; specifies continuation strategy
%                            to be used;
% param.continuation_root,   can be used to specify the continuation root, see
%                            [3, Sec. 2] for building the rational Krylov
%                            space. If not set, used defaults as listed in [2].
% param.continuation_solve,  function handle of the form (default) 
%                            @(AB, nu, mu, x, param) AB.solve(nu, mu, x);
% param.continuation_m,      4 (default); integer which may be used 
%                            by param.continuation_solve; 
% param.continuation_bounds, 0 (default), 1; whether to compute the norms of
%                            \widehat\vf_{j+1} and \widehat\ve_{j+1} used to
%                            produce the bounds related to [3, eq. (3.18)]
%                            and [3, eq. (3.20)]. Works only if param.p = 1,
%                            and param.continuation = 'near-optimal';
%
% The following param field is used to simulate the parallel execution with p
% parallel proces, in the sense that there is less freedom for choosing the
% continuation vectors; cf. [3, Figure 4.2]. It can thus be used to assist
% further research, for instance, when designing approximate linear system
% solvers for the near-optimal continuation strategy for a particular
% application.
% 
% param.p,                   1 (default); integer used to simulate number 
%                            of parallel processors.
%
% If a matrix pencil (A, B) is provided by the structure AB,
% type (iii) for calling the software, it should provide the following.
%
% AB.multiply,  function handle which takes as arguments
%               @(rho, eta, x) and returns rho*A*x-eta*B*x;
% AB.solve,     function handle which takes as arguments
%               @(nu, mu, x) and returns (nu*A-mu*B)\x;
% AB.isreal,    true/false flag specifying whether the pencil is
%               real-valued or not.
%
% Output arguments are:
%
%   - V,      N-by-(m+n) matrix spanning the rational Krylov space;
%   - K, H,   (m+n)-by-(m+n-1) upper-Hessenberg matrices.
%   - param,  the structure from above with some additional
%             information, of interest to the developers.
%
% References:
%
% [1] M. Berljafa and S. G{\"u}ttel. Generalized rational Krylov
%     decompositions with an application to rational approximation,
%     SIAM J. Matrix Anal. Appl., 36(2):894--916, 2015.
%
% [2] M. Berljafa and S. G{\"u}ttel. The RKFIT algorithm for
%     nonlinear rational approximation, MIMS EPrint 2015.38,
%     Manchester Institute for Mathematical Sciences, The University of
%     Manchester, UK, 2015.
%
% [3] M. Berljafa and S. G{\"u}ttel. Parallelization of the rational Arnoldi
%     algorithm, MIMS EPrint 2016.32, Manchester Institute for Mathematical
%     Sciences, The University of Manchester, UK, 2016.
%
% [4] A. Ruhe. Rational Krylov: A practical algorithm for large sparse
%     nonsymmetric matrix pencils, SIAM J. Sci. Comput., 19(5):1535--1551,
%     1998.
%
% [5] A. Ruhe. The rational Krylov algorithm for nonsymmetric eigenvalue
%     problems. III: Complex shifts for real matrices, BIT,
%     34:165--176, 1994.

[AB, V, K, H, N, j, m, rerun, param] = parse_argin(varargin{:});

if param.continuation_bounds && param.p == 1
  W    = zeros(N, m+1);
  Fhat = zeros(1, m);
  fhat = zeros(1, m);
end
  
%R = [0*K(:, 1) 0*K]; % SG: doesn't work for xi = zeros(1,0)
R = zeros(size(K,1),size(K,2)+1); 
T = 0*K;
xi  = param.xi;
nu  = param.moebius(1, :);
mu  = param.moebius(2, :);
rho = param.moebius(3, :);
eta = param.moebius(4, :);

realopt = param.real;

if rerun
  rerun_krylov;
else
  run_krylov;
end

param.R = R;
param.T = T;

param.moebius(1, :) = nu;
param.moebius(2, :) = mu;
param.moebius(3, :) = rho;
param.moebius(4, :) = eta;

if param.continuation_bounds && param.p == 1
  W(:, 1)    = V(:, 1);  
  param.W    = W;
  param.Fhat = Fhat;
  param.fhat = fhat;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Nested functions. %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function run_krylov
  bd = false;
  bd_tol = eps(1);

  % If expanding the array of poles is short.
  offset = j-1;

  % Starting vector.
  if j == 1
    R(1, 1) = param.induced_norm(V(:, 1));
    V(:, 1) = V(:, 1)/R(1, 1);
  end

  if param.waitbar, hwb = waitbar(0, 'rat\_krylov'); end
  while j <= m      
    % Computing the continuation combination.
    continuation_pair;

    % Compute new vector.
    % Continuation vector.
    w = V(:, 1:j)*T(1:j, j);
    
    % Polynomial part.
    w = AB.multiply(rho(j-offset), eta(j-offset), w);            
    
    % Rational part, with eventual refinement.
    if param.refinement, ww = w; end
    w = AB.solve(nu(j-offset), mu(j-offset), w);
    if param.refinement
      ww = ww - AB.multiply(nu(j-offset), mu(j-offset), w);
      w = w + AB.solve(nu(j-offset), mu(j-offset), ww);
    end

    if param.continuation_bounds
      W(:, j+1) = w;
    end
    
    % Orthogonalization.
    if isreal(xi(j-offset)) || realopt == 0
      switch param.orth
       case 'CGS'
        hh = param.inner_product(w, V(:, 1:j));
        w = w - V(:, 1:j)*hh;
        H(1:j, j) = hh;
        if param.reorth == 1
          hh = param.inner_product(w, V(:, 1:j));
          w = w - V(:, 1:j)*hh;
          H(1:j, j) = H(1:j, j) + hh;
        end
        % 'CGS' end
       case 'MGS'
        for reo_i = 1:j
          hh(1) = param.inner_product(w, V(:, reo_i));
          w = w - V(:, reo_i)*hh(1);
          H(reo_i, j) = H(reo_i, j) + hh(1);
        end
        if param.reorth == 1
          for reo_i = 1:j
            hh(1) = param.inner_product(w, V(:, reo_i));
            w = w - V(:, reo_i)*hh(1);
            H(reo_i, j) = H(reo_i, j) + hh(1);
          end
        end
        % 'MGS' end
      end % switch param.orth
      H(j+1, j) = param.induced_norm(w);
      V(:, j+1) = w/H(j+1, j);

      if abs(H(j+1, j)) < bd_tol*norm(H(1:j, j)), bd = true; break; end

      % Setting the decomposition.
      R(1:j+1, j+1) = H(1:j+1, j);
      
      K(1:j+1, j) = H(1:j+1, j);
      K(1:j+1, j) = nu(j-offset)*K(1:j+1, j) - ...
          rho(j-offset)*[T(1:j, j); 0];
      H(1:j+1, j) = mu(j-offset)*H(1:j+1, j) - ...
          eta(j-offset)*[T(1:j, j); 0];

    else
      V(:, j+1) = real(w);
      V(:, j+2) = imag(w);

      switch param.orth % for the real part
       case 'CGS'
        hh = param.inner_product(V(:, j+1), V(:, 1:j));
        V(:, j+1) = V(:, j+1) - V(:, 1:j)*hh;
        H(1:j, j) = hh;
        if param.reorth == 1
          hh = param.inner_product(V(:, j+1), V(:, 1:j));
          V(:, j+1) = V(:, j+1) - V(:, 1:j)*hh;
          H(1:j, j) = H(1:j, j) + hh;
        end
        % 'CGS' end
       case 'MGS'
        for reo_i = 1:j
          hh(1) = param.inner_product(V(:, j+1), V(:, reo_i));
          V(:, j+1) = V(:, j+1) - V(:, reo_i)*hh(1);
          H(reo_i, j) = H(reo_i, j) + hh(1);
        end
        if param.reorth == 1
          for reo_i = 1:j
            hh(1) = param.inner_product(V(:, j+1), V(:, reo_i));
            V(:, j+1) = V(:, j+1) - V(:, reo_i)*hh(1);
            H(reo_i, j) = H(reo_i, j) + hh(1);
          end
        end
        % 'MGS' end
      end % switch param.orth % for the real part
      H(j+1, j) = param.induced_norm(V(:, j+1));
      V(:, j+1) = V(:, j+1)/H(j+1, j);

      if abs(H(j+1, j)) < bd_tol*norm(H(1:j, j)), bd=true; break; end

      j = j+1;
      switch param.orth % for the imag part
       case 'CGS'
        hh = param.inner_product(V(:, j+1), V(:, 1:j));
        V(:, j+1) = V(:, j+1) - V(:, 1:j)*hh;
        H(1:j, j) = hh;
        if param.reorth == 1
          hh = param.inner_product(V(:, j+1), V(:, 1:j));
          V(:, j+1) = V(:, j+1) - V(:, 1:j)*hh;
          H(1:j, j) = H(1:j, j) + hh;
        end
        % 'CGS' end
       case 'MGS'
        for reo_i = 1:j
          hh(1) = param.inner_product(V(:, j+1), V(:, reo_i));
          V(:, j+1) = V(:, j+1) - V(:, reo_i)*hh(1);
          H(reo_i, j) = H(reo_i, j) + hh(1);
        end
        if param.reorth == 1
          for reo_i = 1:j
            hh(1) = param.inner_product(V(:, j+1), V(:, reo_i));
            V(:, j+1) = V(:, j+1) - V(:, reo_i)*hh(1);
            H(reo_i, j) = H(reo_i, j) + hh(1);
          end
        end
        % 'MGS' end
      end % switch param.orth % for the imag part
      H(j+1, j) = param.induced_norm(V(:, j+1));
      V(:, j+1) = V(:, j+1)/H(j+1, j);

      if abs(H(j+1, j)) < bd_tol*norm(H(1:j, j)), bd=true; break; end

      % Setting the decomposition.      
      R(1:j+1, j:j+1) = H(1:j+1, j-1:j);

      rmu = real(mu(j-1-offset)); imu = imag(mu(j-1-offset));
      cmu = [rmu imu; -imu rmu];
      rcnt = [real(T(1:j-1, j-1)) imag(T(1:j-1, j-1)); 0 0; 0 0];           
      
      K(1:j+1, j-1:j) = H(1:j+1, j-1:j);
      K(1:j+1, j-1:j) = K(1:j+1, j-1:j)*nu(j-1-offset) - ...
          rho(j-1-offset)*rcnt;
      H(1:j+1, j-1:j) = H(1:j+1, j-1:j)*cmu - ...
          eta(j-1-offset)*rcnt;
    end % realopt

    if param.waitbar, waitbar(j/m, hwb), end
    j = j+1;
  end % while j <= m
  
  if param.waitbar, close(hwb), end
  
  if bd == true
    warning(['rat_krylov: ''lucky breakdown'' occured at ' ...
             'iteration ' num2str(j)]);
    V = V(:, 1:j);
    K = K(1:j, 1:j-1);
    H = H(1:j, 1:j-1);
  end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Nested functions of rat_krylov. %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function continuation_pair    
  
  j_loc = 1+floor((j-1)/param.p)*param.p;

  if ~isnan(param.continuation_root)
    if isinf(param.continuation_root)
      rho(j-offset) = 0;
      eta(j-offset) = 1;
    else
      rho(j-offset) = 1;
      eta(j-offset) = param.continuation_root;
    end 
  end
    
  switch param.continuation
   case 'almost-last'
    j_loc = max(1, j-param.p+1);
    T(1:j, j) = zeros(j, 1); T(j_loc, j) = 1;
    
   case 'last'
    T(1:j, j) = zeros(j, 1); T(j_loc, j) = 1;
    
   case 'ruhe'
    if j_loc == 1, T(j_loc, j) = 1; else
      [Q, ~] = qr(K(1:j_loc, 1:j_loc-1)*mu(j-offset) - ...
                  H(1:j_loc, 1:j_loc-1)*nu(j-offset));
      T(1:j_loc, j) = Q(:, end);
    end

   case 'near-optimal'
    % Get auxiliary continuation vector. (Other choices may be considered.)
    if j_loc == 1, T(j_loc, j) = 1; else
      [Q, ~] = qr(K(1:j_loc, 1:j_loc-1)*mu(j-offset) - ...
                  H(1:j_loc, 1:j_loc-1)*nu(j-offset));
      T(1:j_loc, j) = Q(:, end);
    end
    
    % Approximate next vector.
    vwhat = AB.multiply(rho(j-offset), eta(j-offset), V(:, 1:j_loc)*T(1:j_loc, j));
    
    if param.continuation_bounds && param.p == 1
      vshat = vwhat;
    end
    
    vwhat = param.continuation_solve(AB, nu(j-offset), mu(j-offset), ...
                                     vwhat, param);    
    
    if param.continuation_bounds && param.p == 1
      % Residual according to [3, eq. (3.8)].
      vshat = AB.multiply(nu(j-offset), mu(j-offset), vwhat) - vshat;
    end
    
    % CGS with reorthogonalization.
    vchat  = param.inner_product(vwhat, V(:, 1:j_loc));
    vvhat  = vwhat - V(:, 1:j_loc)*vchat;
    
    vchat2 = param.inner_product(vvhat, V(:, 1:j_loc));
    vvhat  = vvhat - V(:, 1:j_loc)*vchat2;    
    vchat  = vchat + vchat2;
    
    vchat(end+1, 1) = param.induced_norm(vvhat);
    vvhat = vvhat/vchat(end);
    
    % Form approximate pencil, cf. [3, eq. (3.4)].
    K(1:j_loc+1, j) = nu(j-offset)*vchat - rho(j-offset)*T(1:j_loc+1, j);  
    H(1:j_loc+1, j) = mu(j-offset)*vchat - eta(j-offset)*T(1:j_loc+1, j);

    if param.continuation_bounds && param.p == 1
      % Scaled error given by [3, eq. (3.12)].
      vfhat = -AB.solve(nu(j-offset), mu(j-offset), vshat)/vchat(end);      
      % Corresponding norm, and norm of the local projection.      
      Fhat(j) = abs(param.induced_norm(vfhat));
      fhat(j) = abs(param.induced_norm(...
          [V(:, 1:j) vvhat] * ...
          param.inner_product(vfhat, [V(:, 1:j) vvhat])));
    end    
    
    ind = 1;
    [XX, DD] = eig(H(1:j_loc, [1:j_loc-1 j]), K(1:j_loc, [1:j_loc-1 j]));
  
    if isinf(DD(ind, ind))
      rho(j-offset) = 0;
      eta(j-offset) = 1;
    else
      rho(j-offset) = 1;
      eta(j-offset) = DD(ind, ind);
    end
    
    if isreal(mu(j-offset)) && ...
          isreal(nu(j-offset)) && ...
          isreal(H(1:j_loc, 1:j_loc-1)) && ...
          isreal(K(1:j_loc, 1:j_loc-1)) && ...
          ~isreal(eta(j-offset))
      [~, ind] = min(abs(imag(diag(DD)))); % Not looking for smallest
                                           % imaginary part in relative sense
                                           % to avoid zero.
      eta(j-offset) = DD(ind, ind);
      eta(j-offset) = real(eta(j-offset));
      XX = real(XX);
    end
    
    % Scaling factor given by [3, eq. (3.6)]. (Can actually be excluded.)
    gamma = XX(end, ind)*(rho(j-offset)*H(j_loc+1, j) - ...
                          eta(j-offset)*K(j_loc+1, j));
    
    LL = mu(j-offset)*K(1:j_loc, [1:j_loc-1 j]) - ...
         nu(j-offset)*H(1:j_loc, [1:j_loc-1 j]);
    
    % Near-optimal continuation vector given by [3, eq. (3.7)]. 
    T(1:j_loc, j) = LL*XX(:, ind); 
    T(1:j_loc, j) = T(1:j_loc, j)/gamma;
    K(1:j+1, j) = 0;
    H(1:j+1, j) = 0;        
  end % switch param.continuation  
end % continuation_pair
  
end % run_krylov

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function rerun_krylov
  
  while j <= m      
    % Computing the continuation combination.
    if isreal(xi(j)) || realopt == 0
      T(1:j, j) = K(1:j,j)*mu(j) - H(1:j,j)*nu(j);
      %U(1:j, j) = U(1:j, j)/(nu(j)*eta(j)-mu(j)*rho(j));
    elseif realopt == 1 || realopt == 2      
      [HH, KK, QQ, ZZ] = qz(H(j+1:j+2, j:j+1), K(j+1:j+2, j:j+1));

      HH = H(1:j+2, j:j+1)*ZZ; HH(j+1:j+2, :) = QQ*HH(j+1:j+2, :);
      KK = K(1:j+2, j:j+1)*ZZ; KK(j+1:j+2, :) = QQ*KK(j+1:j+2, :);      
      
      T(1:j, j) = KK(1:j,1)*mu(j) - HH(1:j,1)*nu(j); j = j+1;
      T(1:j, j) = KK(1:j,2)*mu(j) - HH(1:j,2)*nu(j); j = j-1;      
    end
    
    % Next vector and orthogonalization.
    if isreal(xi(j)) || realopt == 0            
      w = V(:, 1:j)*T(1:j, j);    
      w = AB.multiply(rho(j), eta(j), w);
      w = AB.solve(nu(j), mu(j), w);
      
      r = eta(j)*K(1:j+1, j) - rho(j)*H(1:j+1, j);
      %r = r/(nu(j)*eta(j)-mu(j)*rho(j));
      
      % MGS simulation.
      V(:, j+1) = w - V(:, 1:j)*r(1:j, 1);
      V(:, j+1) = V(:, j+1)/r(j+1, 1);     
    else
      w = V(:, 1:j)*T(1:j, j);
      w = AB.multiply(rho(j), eta(j), w);
      w = AB.solve(nu(j), mu(j), w);      

      r = eta(j)*KK(1:j+1, 1) - rho(j)*HH(1:j+1, 1);
      
      % MGS simulation. 
      V(:, j+1) = w - V(:, 1:j)*r(1:j, 1);
      V(:, j+1) = V(:, j+1)/r(j+1, 1);
      
      j = j+1;
      
      w = V(:, 1:j)*T(1:j, j);
      w = AB.multiply(rho(j), eta(j), w);
      w = AB.solve(nu(j), mu(j), w);
      
      r = eta(j)*KK(1:j+1, 2) - rho(j)*HH(1:j+1, 2);
      
      % MGS simulation. 
      V(:, j+1) = w - V(:, 1:j)*r(1:j, 1);
      V(:, j+1) = V(:, j+1)/r(j+1, 1);
            
      V(:, j:j+1) = V(:, j:j+1)*QQ;
      if realopt == 1, V(:, j:j+1) = real(V(:, j:j+1)); end
    end
    
    j = j+1;
  end % while j <= m
end % rerun_krylov

end % rat_krylov



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [AB, V, K, H, N, j, m, ...
          rerun, param] = parse_argin(varargin)
%PARSE_ARGIN    Process the input argument list to rat_krylov.

  msg = 'rat_krylov:parse_argin: ';

  assert(nargin >= 3, [msg 'more parameters needed']);

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%% Defines AB, N, and j. %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  if isstruct(varargin{1})
    AB = varargin{1};
    V_id = 2;
  else
    if size(varargin{2}, 1) == size(varargin{2}, 2) && ...
          size(varargin{2}, 2) ~= 1
      B_id = 2;
      V_id = 3;
    else
      B_id = 0;
      V_id = 2;
    end
  end

  [N, j] = size(varargin{V_id});

  if ~isstruct(varargin{1})
    assert(size(varargin{1}, 1) == N, [msg 'A has to be N-by-N'])
    assert(size(varargin{1}, 2) == N, [msg 'A has to be N-by-N'])
    if B_id == 0
      AB.multiply = @(rho, eta, x) rho*(varargin{1}*x) - eta*x;
      AB.solve    = @(nu, mu, x) (nu*varargin{1}-mu*speye(N))\x;
      AB.isreal = isreal(varargin{1});
    else
      assert(size(varargin{2}, 1) == N, [msg 'B has to be N-by-N'])
      assert(size(varargin{2}, 2) == N, [msg 'B has to be N-by-N'])
      AB.multiply = @(rho, eta, x) ...
          rho*(varargin{1}*x) - eta*(varargin{2}*x);
      AB.solve    = @(nu, mu, x) ...
          (nu*varargin{1}-mu*varargin{2})\x;
      AB.isreal = isreal(varargin{1}) && isreal(varargin{2});
    end
  end

  if ~isfield(AB, 'isreal'), AB.isreal = false; end

  if ~all(isfield(AB, {'multiply', 'solve'}))
    assert(any(isfield(AB, {'A', 'B'})), ...
           [msg 'the pencil is null'])
    if ~isfield(AB, {'B'})
      assert(size(AB.A, 1) == N, [msg 'A has to be N-by-N'])
      assert(size(AB.A, 2) == N, [msg 'A has to be N-by-N'])
      AB.isreal = isreal(AB.A);
      AB.multiply = @(rho, eta, x) rho*(AB.A*x) - eta*x;
      AB.solve    = @(nu, mu, x) (nu*AB.A-mu*speye(N))\x;
    elseif ~isfield(AB, {'A'})
      assert(size(AB.B, 1) == N, [msg 'B has to be N-by-N'])
      assert(size(AB.B, 2) == N, [msg 'B has to be N-by-N'])
      AB.isreal = isreal(AB.B);
      AB.multiply = @(rho, eta, x) rho*x - eta*(AB.B*x);
      AB.solve    = @(nu, mu, x) (nu*speye(N)-mu*AB.B)\x;
    else
      assert(size(AB.A, 1) == N, [msg 'A has to be N-by-N'])
      assert(size(AB.A, 2) == N, [msg 'A has to be N-by-N'])
      assert(size(AB.B, 1) == N, [msg 'B has to be N-by-N'])
      assert(size(AB.B, 2) == N, [msg 'B has to be N-by-N'])
      AB.isreal = isreal(AB.A) && isreal(AB.B);
      AB.multiply = @(rho, eta, x) rho*(AB.A*x) - eta*(AB.B*x);
      AB.solve    = @(nu, mu, x)(nu*AB.A-mu*AB.B)\x;
    end
  end
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%% Defines V, K, H, xi, m, and rerun. %%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  if isstruct(varargin{end}),   param_given = 1;
  elseif ischar(varargin{end}), param_given = 1;
  else                        param_given = 0; end

  rat_krylov_mode = nargin-V_id-param_given;
  rerun = 0;

  switch rat_krylov_mode
   case {1} % Rational Krylov.
    xi = varargin{V_id+1};
    m  = length(xi);
    assert(j==1, [msg 'b has to be N-by-1'])
    V = zeros(N, m+1);
    V(:, 1) = varargin{V_id};
    K = zeros(m+1, m);
    H = zeros(m+1, m);

   case {2} % Rerunning.
    rerun = 1;
    K = varargin{V_id+1};
    H = varargin{V_id+2};
    assert(j==1, [msg 'b has to be N-by-1'])
    assert(size(K, 1) == size(H, 1), ...
           [msg 'K and H need to be of equal size'])
    assert(size(K, 2) == size(H, 2), ...
           [msg 'K and H need to be of equal size'])
    assert(size(K, 1) == size(K, 2)+1, ...
           [msg 'K has to be of size (m+1)-by-m'])
    m = size(K, 2);
    V = zeros(N, m+1);

    % Allow for multiple precision rerunning.
    if isa(K, 'sym') || isa(H, 'sym') || isa(varargin{V_id}, 'sym')
      V = sym(V);   
    elseif isa(K, 'mp') || isa(H, 'mp') || isa(varargin{V_id}, 'mp')
      V = mp(V);
    end

    V(:, 1) = varargin{V_id};
    xi = util_pencil_poles(K, H);

   case {3} % Expanding.
    xi = varargin{V_id+3};
    m  = length(xi);
    assert(size(varargin{V_id+1}, 1) == j, ...
           [msg 'K has to be of size j-by-(j-1)'])
    assert(size(varargin{V_id+1}, 2) == j-1, ...
           [msg 'K has to be of size j-by-(j-1)'])
    assert(size(varargin{V_id+2}, 1) == j, ...
           [msg 'K has to be of size j-by-(j-1)'])
    assert(size(varargin{V_id+2}, 2) == j-1, ...
           [msg 'K has to be of size j-by-(j-1)'])
    V = zeros(N, m+j);
    K = zeros(m+j, m+j-1);
    H = zeros(m+j, m+j-1);
    V(:, 1:j) = varargin{V_id};
    K(1:j, 1:j-1) = varargin{V_id+1};
    H(1:j, 1:j-1) = varargin{V_id+2};
    m = m+j-1;

   otherwise
    error([msg 'invalid number of input arguments']);
  end
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%% Defines param. %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  if ~param_given || ischar(varargin{end})
    param = struct;
    if strcmp('real', varargin{end}), param.real = 1; end
  else
    param = varargin{end};
  end

  if isfield(param, 'orth')
    if ~strcmp(param.orth, 'CGS') && ~strcmp(param.orth, 'MGS')
      warning('Warning: param.orth invalid. Setting to default.')
      param.orth = 'MGS';
    end
  else, param.orth = 'MGS'; end

  if isfield(param, 'reorth')
    if param.reorth ~= 0 && param.reorth ~= 1
      warning('Warning: param.reorth invalid. Setting to default.')
      param.reorth = 1;
    end
  else, param.reorth = 1; end

  if isfield(param, 'real')
    if param.real ~= 1 && param.real ~= 0
      warning('Warning: param.real invalid. Setting to default.')
      param.real = 0;
    end
  else, param.real = 0; end

  if isfield(param, 'refinement')
    if param.refinement ~= 0 && param.refinement ~= 1
      warning('Warning: param.refinement invalid. Setting to default.')
      param.refinement = 0;
    end
  else, param.refinement = 0; end


  if rat_krylov_mode == 2 && size(H, 2) > 1 && any(diag(H, -2))
    param.real = 1;
  end

  if param.real && (~AB.isreal || ~isreal(V(:, 1)))
    if rat_krylov_mode == 2, param.real = 2;
    else                     param.real = 0;
      warning(['Warning: Ignoring ''real'' option with complex' ...
               ' matrices/vectors.']);
    end
  end
  if param.real && ~canonical_cplx(xi)
    param.real = 0;
    warning(['Warning: Ignoring ''real'' option as conj(poles(j))' ...
             ' == poles(j+1) failed']);
  end

  param.xi = xi;
  param.moebius = poles_to_moebius(xi);

  if ~isfield(param, 'inner_product')
    param.inner_product = @(x, y) y'*x;
  end
  param.induced_norm  = @(x) sqrt(param.inner_product(x, x));
  
  if isfield(param, 'waitbar')
    if param.waitbar ~= 1 && param.waitbar ~= 0
      warning('Warning: param.waitbar invalid. Setting to default.')
      param.waitbar = 0;
    end
  else, param.waitbar = 0; end    
  
  if ~isfield(param, 'continuation')
    param.continuation = 'ruhe';
  end

  if ~isfield(param, 'continuation_root')
    param.continuation_root = NaN;
  end  
  
  if ~isfield(param, 'continuation_solve')
    param.continuation_solve = @(AB, nu, mu, x, param) AB.solve(nu, mu, x);
  end  
  
  if ~isfield(param, 'continuation_bounds')
    param.continuation_bounds = 0;
  end

  if ~isfield(param, 'continuation_m')
    param.continuation_m = 4;
  end  

  % To simulate parallel execution. 
  
  if ~isfield(param, 'p')
    param.p = 1;
  end  

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    

  assert(rerun || m < N, [msg 'number of poles cannot be greater' ...
                    ' or equal to the size of the problem']);
  assert(size(xi, 1) == 1 , ...
         [msg 'xi must be a 1-by-m row vector']);
  % assert(all(~any(tril(K(3:end, 1:end-1)))), ...
  %         [msg 'H is not upper-Hessenberg']);
  % assert(all(~any(tril(H(4:end, 1:end-1)))), ...
  %         [msg 'K is not quasi-upper-Hessenberg']);

end % parse_argin



function y = canonical_cplx(xi)
% CANONICAL_CPLX Check if the poles xi are ordered canonically.

  m = length(xi);
  y = 1;
  j = 1;
  while j <= m
    if isreal(xi(j)) || isinf(double(xi(j)))
      j = j+1;
    else
      if j == m || (j < m && xi(j+1) ~= conj(xi(j)))
        y = 0;
      end % if
      j = j+2;
    end % if
  end % while

end


function moebius = poles_to_moebius(xi)
% POLES_TO_MOEBIUS Moebius transformation with poles xi.
%
% Finite xi is replaced with (nu, mu) := (1, xi) and
% (rho, eta) := (0, 1),  and xi = inf is replaced by
% (nu, mu) := (0, 1) and (rho, eta) := (1, 0).

  nu  = ones(1, length(xi));
  mu  = xi;
  rho = zeros(1, length(xi));
  eta = ones(1, length(xi));

  %nu(isinf(xi))  = 0; % does not work with vpa
  nu(abs(xi)==inf) = 0;
  %mu(isinf(xi))  = 1; % does not work with vpa
  mu(abs(xi)==inf) = 1;
  rho(abs(xi)>1) = 1;
  eta(abs(xi)>1) = 0;
  
  moebius = [nu; mu; rho; eta];

end
