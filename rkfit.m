function [poles,ratfun,misfit,K,H,V,coeffs] = rkfit(F,A,v,poles,realflag)
%RKFIT    Rational least squares fitting via rational Krylov.
%
% Attempts to find an improved set of poles for a rational function r(z) 
% such that || F*v - r(A)*v ||_2 is small and, after iterating this
% function, hopefully minimal among all rational functions of type (m,m).
%
% Calling syntax: [poles,ratfun,misfit,K,H,V,coeffs] = rkfit(F,A,v,poles,realflag)
%
% The inputs are as follows:
%
%   - F is either a N-by-N matrix or a cell-array of N-by-N matrices;
%   - A is a square N-by-N matrix;
%   - v is an N-by-1 column vector;
%   - poles is a 1-by-m row vector of initial poles;
%   - realflag (optional). If set to 'real' the computed poles will appear
%     in perfectly complex conjugate pairs. This is useful when F, A, and v
%     are real and the initial poles appear in complex conjugate pairs
%     (use CPLXSORT to order the poles accordingly).
%
% The outputs are as follows:
%
%   - poles is a 1-by-m row vector of the new (improved) poles;
%   - ratfun is a function handle with variable length input argument list:
%     * if ratfun is called with one input argument, ratfun(Z),
%       then the scalar rational approximant will be evaluated for each
%       entry in Z,
%     * if ratfun is called with input arguments, ratfun(A,v),
%       then r(A)*v is evaluated as a matrix function times a vector 
%       or matrix v;
%   - misfit, which is a real number equal to || F*v - r(A)*v ||_2;
%   - K and H are (m+1)-by-m unreduced upper Hessenberg matrices;
%   - V is an N-by-(m+1) orthonormal rational Krylov basis.
%
% This algorithm is described in
%
% [1] M. Berljafa and S. G\"{u}ttel. Generalized rational Krylov
%     decompositions with an application to rational approximation,
%     MIMS EPrint 2014.59, Manchester Institute for Mathematical
%     Sciences, The University of Manchester, UK, 2014. 

if nargin == 5 && strcmp(realflag,'real'),
    realflag = 'real';
else
    realflag = 'complex';
end

if ~iscell(F)
  F = {F};
end
block = length(F);

if iscell(v) || size(v, 2) ~= 1
  if iscell(v)
    w = zeros(size(A, 1), length(v));
    for j = 1:block    
      w(:, j) = v{j};
    end
  else
    w = v;    
  end 
  v = ones(size(A, 1), 1);
else
  w = ones(size(A, 1), block);
end

poles =  cplxsort(poles(:).');
[V, K, H] = rat_krylov(A, v, poles, realflag);

MM = [];
for j = 1:block
  M = F{j}*V;
  M = V*(V'*M)-M;
  MM = [MM; diag(w(:, j))*M];
end

[~, ~, W] = svd(MM,0); 
W = W(:, end:-1:1);
K = W'*K; H = W'*H;

poles = cplxsort((eig(H(2:end,:), K(2:end,:))));
poles = poles(:).';

% compute coeffs of new rational approximant
if nargout > 1
    [V, K, H] = rat_krylov(A, v, poles, realflag);
    
    coeffs = cell(1, block);
    ratfun = cell(1, block);
    misfit = 0.0;
    
    for j = 1:block
      Fv = F{j}*v;
      coeffs{j} = V'*Fv;
      misfit = misfit + norm(V*coeffs{j} - Fv).^2;
      coeffs{j} = coeffs{j}/norm(v);
      ratfun{j} = @(varargin) rat_eval(K, H, coeffs{j}, varargin);
    end
    
    misfit = sqrt(misfit);
    
    if block == 1
      ratfun = ratfun{1};
      coeffs = coeffs{j};
    end
end
end



function w = rat_eval(K,H,coeffs,varargin)
%RAT_EVAL   Evaluate a rational function expansion associated with K and H.
%
% Calling syntax: w = rat_eval(K,H,coeffs,A)
%                 w = rat_eval(K,H,coeffs,A,v)
%
% The rational function r(A) can either be evaluated pointwise (if only A
% is provided), or as a matrix function r(A)*v if A and v are provided.
%
% This function is used mainly by RKFIT which is described in
%
% [1] M. Berljafa and S. G\"{u}ttel. Generalized rational Krylov
%     decompositions with an application to rational approximation,
%     Manchester EPrint available, 2014.

A = varargin{1}{1};

if length(varargin{1}) == 1  % point wise evaluation
    [mA,nA] = size(A);
    A = spdiags(A(:),0,mA*nA,mA*nA);
    v = ones(length(A),1);
else                         % matrix function times vector
    v = varargin{1}{2};
end

w = zeros(size(v));
for j = 1:size(v,2),
    V = rat_krylov(A, v(:,j), K, H);
    w(:,j) = V*coeffs;
end

if length(varargin{1}) == 1,
    w = reshape(w,mA,nA);
end
end
