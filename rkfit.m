function [xi, ratfun, misfit, out] = rkfit(F, A, b, xi, varargin)
% RKFIT    Rational Krylov fitting.
%
% This functions implements the rkfit algorithm [1, 2] used for
% nonlinear rational least squares fitting.
%
% In the simple case when F and A are N-by-N matrices and b is a
% N-by-1 vector rkfit attempts to find a rational function r such
% that \|F*b-r(A)*b\|_2/\|F*b\|_2 is minimal, using xi as the
% initial guess for the poles of r, and iteratively relocating it.
% See eq. (1.5) in [2] for a precise formulation of the general
% problem.
%
% More general constructions are allowed; F may be a cell array
% specifying a sequence of ell matrices of size N-by-N, and b may
% be a N-by-n block of vectors. Additional weighting matrices can
% be passed, as well as additional parameters as detailed below.
%
% We now list the possible calls. The output is always of the form
% [xi, ratfun, misfit, out], and is therefore excluded from the
% following list.
%
%   rkfit(F, A, b, xi)
%   rkfit(F, A, b, xi, maxit)
%   rkfit(F, A, b, xi, maxit, tol)
%   rkfit(F, A, b, xi, maxit, tol, 'real')
%   rkfit(F, A, b, xi, param)
%
% The input arguments are:
%   - F,      N-by-N matrix or a function handle allowing one to
%             perform matrix matrix multiplication, e.q.,
%             F = @(X) Fmat*X, alternatively, F may be a cell array
%             containing ell of such objects;
%   - A,      N-by-N matrix, or a struct representing a N-by-N
%             pencil, see rat_krylov for details;
%   - b,      N-by-n matrix;
%   - xi,     1-by-m row-vector of starting poles, with m < N;
%   - maxit,  equivalent to param.maxit, see below;
%   - tol,    equivalent to param.tol, see below;
%   - 'real', equivalent to setting param.real = 1, see below;
%   - param,  struct, see below.
%
% The fields of the param structure are listed below. The names are
% self-explanatory, and if not provided, the default is being used.
%
% param.k          0 (default); integer used for the specification
%                  of the type (m+k, m) for the rational
%                  approximants;
% param.D,         of the same data-type as F, used for weighted
%                  least squares;
% param.maxit,     10 (default); natural number specifying the
%                  maximal number of iterations rkfit is allowed to
%                  perform;
% param.real,      0 (default) or 1 to enforce real computation
%                  only if possible;
% param.reduction, 1 (default) or 0; reduction occurs only if the
%                  variable is set to 1, provided that relative
%                  misfit below param.tol has been reached;
% param.return,    'best' (default) or 'last'; which of the rkfuns
%                  to return, judged by the misfit; if reductions
%                  takes places and 'best' is on it will return the
%                  best approximation with lower denominator degree
%                  only if it provides misfit below param.tol,
%                  otherwise the best with higher degree;
% param.safe,      safety parameter for the reduction defaults to
%                  0.1, see eq. (4.3) in [2] for details;
% param.stable,    0 (default) or 1; whether to enforce stability
%                  of poles in the relocation or not;
% param.tol,       tolerance to achieve with the relative misfit,
%                  by default 1e-15.
%
% Other fields for param may be those used by RAT_KRYLOV.
%
% Output arguments are:
%
%   - xi,     1-by-m row-vector containing last obtained poles;
%   - ratfun, rkfun object representing the obtained approximant,
%             or if ell > 1, a cell array of rkfuns;
%   - misfit, row-vector of relative misfits per iteration;
%   - out,    struct containing
%                  tol: param.tol from above for reference;
%       misfit_initial: relative misfit from the starting guess;
%            m_initial: initially chosen m parameter;
%            m_reduced: reduced degree, may be less or equal to m
%            k_initial: initially chosen k parameter;
%            k_reduced: ell-vector of reduced k parameters;
%           xi_initial: vector with initial poles;
%         xi_unreduced: poles after convergence;
%           xi_reduced: roots of approximate GCD;
%         reduction_it: NaN or iteration when reduction occurred;
%            return_it: iteration from which approximant is returned;
%      singular_values: of S during the reduction procedure;
%                    F: F, A and b are provided for reference only
%                    A: and their representation may be
%                    b: different from the inputted one;
%                    V: orthonormal basis from last iteration.
%
% Example: Fitting a rational function (A + 4I)\(A - 2.5I) with pole at -4
%          starting with an initial guess of m=2 infinite poles and k=2.
%
% A = diag(1:50); I = eye(50); b = ones(50, 1);
% [xi, ratfun, misfit, out] = rkfit((A+4*I)\(A-2.5*I), A, b, [inf, inf], ...
%                                   struct('real',      1, ...
%                                          'tol',       1e-10, ...
%                                          'maxit',     10, ...
%                                          'k',         2));
%
% >> xi
% xi =
%    -4.0000
% >> ratfun
% ratfun =
%     RKFUN object of type (1,1).
%     Real Hessenberg pencil (H,K) of size 2-by-1.
%     coeffs = [ 0.676, 0.258 ]
%
% References
%
% [1] M. Berljafa and S. G{\"u}ttel. Generalized rational Krylov
%     decompositions with an application to rational approximation,
%     SIAM J. Matrix Anal. Appl., 36(2):894--916, 2015.
%
% [2] M. Berljafa and S. G{\"u}ttel. The RKFIT algorithm for
%     nonlinear rational approximation, MIMS EPrint 2015.38,
%     Manchester Institute for Mathematical Sciences,
%     The University of Manchester, UK, 2015.
%
% See also RKFUN, RAT_KRYLOV.

% Undocumented stuff:
% param.relsvd, 0 (default) or 1; 0 means not used

  [F, ...    % Cell array of function handles to fit.
   A, ...    % Square matrix of dimension N*n.
   b, ...    % Starting vector of size N*n-by-1.
   N, ...    % Natural number.
   n, ...    % Natural number such that n <= N.
   ell, ...  % Number of elements in the cell arrays F, D.
   xi, ...   % Vector of size 1-by-m with initial poles.
   D, ...    % Cell array of function handles for weights.
   param ... % Additional parameters.
  ] = parse_argin_rkfit(F, A, b, xi, varargin{:});

  Nn = N*n;
  m = length(xi);

  % Initialize the out structure.
  fill_out;
  
  if nargout == 4
    out.F = F;
    out.A = A;
    out.b = b;
  end
  
  if strcmp('best', param.return), return_best = 1;
  else                             return_best = 0; end  
  
  if return_best
    best_m = m;
    best_misfit = inf;
  end
  
  coeffs = cell(1, ell);
  ratfun = cell(1, ell);

  if param.reduction == 0, reduction_done  = true;
  else                     reduction_done  = false; end  
  % The following flag is used to indicate that the attempted
  % reduction of m eventually reveals there is no need for
  % reduction, in which case we no longer need to iterate.  
  reduction_break = false;
  % The following flag is used to compute F*v or D*F*v only once,
  % at the very first iteration.
  target_computed = false;
  poles_perturbation = [];

  % The loop does at most param.maxit+1 iterations, which means at  
  % most param.maxit relocations of poles. The eventual
  % (param.maxit+1)st iteration only compues the final rational
  % Krylov space and projects the (D*)F*v onto it.
  it = 1;
  while it <= param.maxit+1
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Construct the search and target spaces. %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % See Section 3.1 in [2].
    if param.k > 0, xi = [xi inf*ones(1, param.k)]; end
    [V, K, H] = rat_krylov(A, b, xi, param);
    if param.k < 0, [~, ~, Qnum] = util_hh2th(K, H); end
    if param.k >= 0, Qnum = eye(m+1+param.k);
    else             Qnum = Qnum(1:m+1+param.k, :)'; end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Compute the approximant. %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Compute F*v only once, and also the norm.
    if ~target_computed
      Fv = repmat(0, Nn, ell);
      nrm_Fv = zeros(1, ell);
      for j = 1:ell
        Fv(:, j) = F{j}(V(:, 1));
        if param.weight, Fv(:, j)  = D{j}(Fv(:, j)); end
        nrm_Fv(j) = norm(Fv(:, j), 'fro');
      end

      nrm_denom = norm(nrm_Fv, 'fro');
      target_computed = true;
    end

    % Project F*v onto the target space, and compute the current
    % misfit according to eq. (1.5) in [2].
    local_misfit = 0;
    for j = 1:ell
      if ~param.weight
        coeffs{j} = Qnum*(Qnum'*(V'*Fv(:, j)));
        Ev = V*coeffs{j}-Fv(:, j);
      else
        DV = D{j}(V);
        coeffs{j} = Qnum*(Qnum'*((DV'*DV)\(DV'*Fv(:, j))));
        Ev = DV*coeffs{j}-Fv(:, j);
      end
      ratfun{j} = rkfun(K, H, coeffs{j}, param.k);
      local_misfit = local_misfit + norm(Ev, 'fro')^2;
    end
    local_misfit = sqrt(local_misfit)/nrm_denom;
    
    misfit(it) = local_misfit;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Save if best solution thus far. %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if return_best 
      save_rf = or(local_misfit < best_misfit, ...
                   and(m < best_m, local_misfit < param.tol));
      if save_rf
        best_m  = m;
        best_ratfun = ratfun;
        best_coeffs = coeffs;
        best_misfit = local_misfit;
        best_V = V; best_K = K; best_H = H; best_xi = xi;     
        out.return_it = it;
      end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Check misfit and maybe stop. %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if it == param.maxit+1 || m == 0, break, end
    if local_misfit < param.tol && reduction_done, break, end

    %%%%%%%%%%%%%%%%%%%%%%%%%
    % Build the SVD matrix. %
    %%%%%%%%%%%%%%%%%%%%%%%%%

    % Compute matrix S for the SVD.
    if it == 1, S = repmat(0, Nn*ell, m+1); end
    if it > 1 && size(S, 2) > m+1, S = S(:, 1:m+1); end
    for j = 1:ell
      jb = (j-1)*Nn + 1;
      je = j*Nn;
      % Compute F*V_search.
      S(jb:je, :) = F{j}(V(:, 1:m+1));
    end

    if param.relsvd, [UU, SS, WW] = svd(S,0); S = UU; end

    for j = 1:ell
      jb = (j-1)*Nn + 1;
      je = j*Nn;
      % Compute (I-V_target*V_target')*F*V_search.
      S(jb:je, :) = S(jb:je, :) - (V*Qnum)*((Qnum'*V')*S(jb:je, :));
      % Eventually apply weighting.
      if param.weight, S(jb:je, :) = D{j}(S(jb:je, :)); end
    end

    if local_misfit < param.tol
      [~, s, W] = svd(S, 0);
      s = diag(s);
      reduce_denominator;
    else
      [~, ~, W] = svd(S, 0);
      W = W(:, end:-1:1);

      if param.relsvd, W = WW*pinv(SS)*W(:, 1); [W, ~] = qr(W); end

      K = W'*K(1:m+1, 1:m);
      H = W'*H(1:m+1, 1:m);

      xi = util_cplxpair((eig(H(2:m+1,1:m), K(2:m+1,1:m))));
      xi = xi(:).';
      if param.stable
        xi(real(xi)>0) = ...
            -real(xi(real(xi)>0)) + ...
            1i*imag(xi(real(xi)>0));
      end

    end % while it <= param.maxit+1

    if reduction_break, break, end

    it = it + 1;
  end

  if return_best && out.return_it ~= it
    coeffs = best_coeffs;
    ratfun = best_ratfun;
    V = best_V; K = best_K; H = best_H; xi = best_xi;    
    clear best_V best_K best_H best_xi;
  else
    out.return_it = it;
  end
  
  xi = xi(1:m);

  if param.reduction ~= 0, reduce_numerator; end

  if ell == 1, coeffs = coeffs{1}; ratfun = ratfun{1}; end

  out.misfit_initial = misfit(1);
  if length(misfit) >= 2, misfit = misfit(2:end);
  else                    misfit = [];            end

  if nargout == 4, out.V = V; end    
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Nested functions. %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function reduce_denominator

% Ensure reduction can be performed only once.
  reduction_done = true;
  out.singular_values = s;

  % Determine m_star.
  % See eq. (4.3) in [2].
  m_star = length(s) - sum(s < param.tol*param.safe*nrm_denom);

  % Reduction may not be needed if m was chosen ``correctly''.
  if m_star >= m, reduction_break = true; return, end

  if param.k + m_star < 0, m_star = -param.k; end
  defect = m-m_star+1;

  W = W(:, end:-1:1);
  if param.relsvd, W = WW*pinv(SS)*W; [W, ~] = qr(W); end

  KK = W'*K(1:m+1, 1:m);
  HH = W'*H(1:m+1, 1:m);

  % The following is not necessary. It may enhance the performance
  % of the relocation of poles for the case that F is not a
  % rational function, c.f. [Lemma 4.3, 2].
  if param.reduction == 1
    [~, ~, WW] = svd([HH(1+defect:m+1,:); KK(1+defect:m+1,:)]);    
    WW = WW(:, end:-1:1);
    KK = KK*WW;
    HH = HH*WW;
  end

  xi_new = eig(HH(1+defect:m+1,defect:m), KK(1+defect:m+1,defect:m));
  if param.stable
    xi_new(real(xi_new)>0) = ...
        -real(xi_new(real(xi_new)>0)) + ...
        1i*imag(xi_new(real(xi_new)>0));    
  end

  xi_new = util_cplxpair(xi_new);

  out.m_reduced = m_star;
  out.xi_unreduced = xi(:).';
  out.xi_reduced   = xi_new(:).';

  out.reduction_it = it;
  xi = xi_new(:).';
  m  = m_star;
end % reduce_denominator


function reduce_numerator

  local_misfit = 0;
  [~, ~, Q] = util_hh2th(K, H);
  for j = 1:ell
    if ~param.weight
      tmp = Q(1:m+param.k+1, :)*coeffs{j};
      nrm_cml = sqrt(cumsum(tmp(end:-1:1).^2));
      % See eq. (4.2) in [2].
      delta_k = sum(nrm_cml <= (param.tol*nrm_denom - ...
                                nrm_denom*misfit(end))/ell);
      out.k_reduced(j) = param.k-delta_k;
      tmp(end-delta_k+1:end) = 0;
      coeffs{j} = Q(1:m+param.k+1, :)'*tmp;
      if param.k-delta_k <= 0
        ratfun{j} = rkfun(K(1:m+1, 1:m), ...
                          H(1:m+1, 1:m), ...
                          coeffs{j}(1:m+1), param.k-delta_k);
      else
        coeffs{j} = coeffs{j}(1:end-delta_k);
        ratfun{j} = rkfun(K(1:end-delta_k, 1:end-delta_k), ...
                          H(1:end-delta_k, 1:end-delta_k), ...
                          coeffs{j}, param.k-delta_k);
      end
      Ev = V(:, 1:length(coeffs{j}))*coeffs{j}-Fv(:, j);
    else
      warning('reduce_numerator not supported with weighting');
      DV = D{j}(V);
      Ev = DV*coeffs{j}-Fv(:, j);
    end
    local_misfit = local_misfit + norm(Ev, 'fro')^2;
  end % j = 1:ell

  if return_best
    misfit(out.return_it) = sqrt(local_misfit)/nrm_denom;
  else 
    misfit(it) = sqrt(local_misfit)/nrm_denom;
  end
end % reduce_numerator

function fill_out
  out.tol = param.tol;
  out.misfit_initial = inf;
  out.m_initial = m;
  out.m_reduced = m;
  out.k_initial = param.k;
  out.k_reduced = param.k*ones(1, ell);
  out.xi_initial   = xi;
  out.xi_unreduced = [];
  out.xi_reduced   = [];
  out.reduction_it = NaN;
  out.return_it    = NaN;
  out.singular_values = [];
end

end % rkfit



function [F, A, b, N, n, ell, xi, D, param] = ...
      parse_argin_rkfit(F, A, b, xi, varargin)
%PARSE_ARGIN_RKFIT    Process the input argument list to rkfit.

  msg = 'parse_argin_rkfit: ';

  param = struct;
  D  = [];

  % Set param.k, param.weight, param.real, and D.
  switch nargin
   case 4+0
    param.k = 0;
    param.weight = false;
    param.real   = 0;

   case 4+1
    if isstruct(varargin{1})
      param = varargin{1};
      if ~isfield(param, 'k'), param.k = 0; end
      if  isfield(param, 'D'), param.weight = true; D = param.D;
      else                     param.weight = false; end
      if ~isfield(param, 'real'), param.real = 0; end
    else
      param.k = 0;
      param.real = 0;
      param.weight = false;
      if ischar(varargin{1})
        if strcmp('real', varargin{1}), param.real = 1; end
      else
        param.maxit = varargin{1};
      end
    end

   case 4+2
    param.k = 0;
    param.real = 0;
    param.weight = false;
    char_id = ischar(varargin{1})+2*ischar(varargin{2});
    if char_id
      if strcmp('real', varargin{char_id}), param.real = 1; end
      param.maxit = varargin{3-char_id};
    else
      param.tol = varargin{2};
      param.maxit = varargin{1};
    end

   case 4+3
    param.k = 0;
    param.real = 0;
    param.weight = false;
    char_id = ischar(varargin{1})+2*ischar(varargin{2})+3* ...
              ischar(varargin{3});
    if strcmp('real', varargin{char_id}), param.real = 1; end
    other_id = [1:char_id-1 char_id+1:3];
    param.tol = varargin{other_id(2)};
    param.maxit = varargin{other_id(1)};

   otherwise
    error([msg 'invalid number of input arguments']);
  end

  N = size(b, 1);
  n = size(b, 2);
  b = b(:);

  assert(n <= N, [msg 'invalid number of starting vectors']);

  if ~iscell(A), A = {A}; end
  if ~iscell(F), F = {F}; end
  if ~iscell(D), D = {D}; end

  len_d = length(D);
  len_f = length(F);

  A = construct_pencil(A, N, n);

  for j = 1:len_f
    if ~isa(F{j}, 'function_handle')
      F{j} = @(X) F{j}*X;
    end
  end

  if param.weight
    for j = 1:len_f
      if ~isa(D{j}, 'function_handle')
        D{j} = @(X) D{j}*X;
      end
    end
  end

  ell = len_f;
  if param.weight, ell = max(ell, len_d); end

  assert(len_d == 1 || len_d == ell, [msg 'invalid ell']);
  assert(len_f == 1 || len_f == ell, [msg 'invalid ell']);

  if ell > 1 && param.weight
    if len_d == 1, for j = 2:ell, D{j} = D{1}; end, end
    if len_f == 1, for j = 2:ell, F{j} = F{1}; end, end
  end

  % So far, F{j} and D{j} represent N-by-N operators, but A
  % has been enlarged to N*n-by-N*n.
  for j = 1:ell
    if n > 1, F{j} = @(X) apply_target(n, F{j}, X); end
    if n > 1, D{j} = @(X) apply_target(n, D{j}, X); end
  end

  xi = util_cplxpair(xi(:).');

  if ~isfield(param, 'maxit'),     param.maxit = 10;       end
  if ~isfield(param, 'tol'),       param.tol = 1e-15;      end
  if ~isfield(param, 'k'),         param.k = 0;            end
  if ~isfield(param, 'relsvd'),    param.relsvd = false;   end
  if ~isfield(param, 'reduction'), param.reduction = 1;    end
  if ~isfield(param, 'safe'),      param.safe = 1e-1;      end
  if ~isfield(param, 'stable'),    param.stable = false;   end
  
  if ~isfield(param, 'return')
    param.return = 'best';  
  else
    valid = or(strcmp('best', param.return), ...
               strcmp('last', param.return));
    if ~valid
      warning(['rkfit: param.return invalid; defaulting to' ...
               ' ''best''']);
      param.return = 'best';
    end
  end  
  
  
end



function AB = construct_pencil(A, N, n)
% Used to construct the pencil structure used by rat_krylov to
% represent the matrix A.

  msg = 'construct_pencil: ';

  C = cell(1, length(A));
  ABisreal = true;

  for i = 1:length(A)
    if isstruct(A{i})
      assert(isfield(A{i}, 'multiply'), ...
             [msg 'A{' num2str(i) '} misses ''mutiply'' field'])
      assert(isfield(A{i}, 'solve'), ...
             [msg 'A{' num2str(i) '} misses ''solve'' field'])
      C{i} = A{i};
      if isfield(A{i}, 'isreal')
        ABisreal = ABisreal && A{i}.isreal;
      else
        ABisreal = false;
      end
    else
      assert(size(A{i}, 1) == N, 'A must be a N-by-N matrix')
      assert(size(A{i}, 2) == N, 'A must be a N-by-N matrix')
      C{i}.multiply = @(rho, eta, x) rho*(A{i}*x) - eta*x;
      C{i}.solve    = @(nu, mu, x) (nu*A{i}-mu*speye(N))\x;
      C{i}.isreal   = isreal(A{i});
      ABisreal = ABisreal && C{i}.isreal;
    end
  end

  if n == 1
    AB = C{1};
  else
    AB.multiply = @(rho, eta, x) ...
        apply_pencil(0, n, C, rho, eta, x);
    AB.solve = @(nu, mu, x) ...
        apply_pencil(1, n, C, nu, mu, x);
    AB.isreal = ABisreal;
  end
end % function construct_pencil



function X = apply_pencil(type, n, C, a, b, X)
% Helper function for construct_pencil used to define multiply and
% solve for the pencil.

  Nn = size(X, 1);
  N = Nn/n;
  assert(mod(Nn, n) == 0, 'cell_multiply: dimension error');
  if n > 1 && length(C) == 1, map1 = @(index) 1;
  else                        map1 = @(index) index; end

  for j = 1:n
    jb = (j-1)*N+1;
    je = j*N;
    if type == 0
      X(jb:je, :) = C{map1(j)}.multiply(a, b, X(jb:je, :));
    else % type == 1
      X(jb:je, :) = C{map1(j)}.solve(a, b, X(jb:je, :));
    end
  end
end % function apply_pencil



function result = apply_target(n, F, X)
% Used to apply F to X exploiting the Kronecker structure of F.

  Nn = size(X, 1);
  N = Nn/n;

  assert(mod(Nn, n) == 0, 'apply_target: dimension error');

  result = repmat(0, size(X, 1), size(X, 2));
  for j = 1:n
    jb = (j-1)*N+1;
    je = j*N;
    result(jb:je, :) = F(X(jb:je, :));
  end
end % function apply_target


% The following function is currently not used.
function index = greedy_match(x, y, m_star, real)
% Used to determine the m_star elements in y closest to some m_star
% elements in x, where x and y are vectors of length n > m_star,
% and index is a {0, 1}-valued vector of length n with m_star
% elements set to 1.

  n = length(x);
  d = abs(repmat(x(:), 1, length(y))-repmat(y(:).', n, 1));
  d = d./(abs(repmat(x(:), 1, length(y)))+abs(repmat(y(:).', n, 1)));
  index = zeros(n, 1);

  while m_star > 0
    [~, c]  = min(min(d));
    [~, r]  = min(d(:, c));
    d(:, c) = inf;
    d(r, :) = inf;
    index(r) = 1;
    if real && imag(x(r)) ~= 0
      r_pair = find(x == x(r)', 1);
      c_pair = find(y == y(c)', 1);
      d(:, c_pair) = inf;
      d(r_pair, :) = inf;
      index(r_pair) = 1;
      m_star = m_star - 2;
    else
      m_star = m_star - 1;
    end
  end

end % function greedy_match
