function [f, f2, fmp] = rat_qr_eval(K, H, c, z)
%RAT_QR_EVAL   Evaluate rational function via QR.
% f = rat_qr_eval(K, H, c, z)
%
% This algorithm is described in
%
% [1] M. Berljafa and S. G\"{u}ttel. Generalized rational Krylov
%     decompositions with an application to rational approximation,
%     MIMS EPrint 2014.59, Manchester Institute for Mathematical
%     Sciences, The University of Manchester, UK, 2014. 

  m = size(H, 2);
  
  f = z;
  for i = 1:size(z,1),
    for j = 1:size(z,2),
      L = z(i,j)*K-H;  
      D = arrayfun(@(col) norm(L(:, col)), 1:m);
      D = floor(log2(D));
      L = L/diag(pow2(D));
      [Q, ~] = qr(L);
      f(i,j) = Q(:, end)'*c/conj(Q(1, end));
    end
  end

  f2 = z;  
  for i = 1:size(z,1)
    for j = 1:size(z,2)      
      L = z(i,j)*K-H;   
      D = arrayfun(@(col) norm(L(:, col)), 1:m);
      D = pow2(floor(log2(D)));
      L = L/diag(D);
      ql = eye(m+1, 1); % last-column
      rr = ql;          % r(k)
      
      for k = 1:m       
        L(1:k+1, k) = L(1:k+1, k) - L(1:k+1, 1:k-1)*(L(1:k+1, 1:k-1)'*L(1:k+1, k));
        L(1:k+1, k) = L(1:k+1, k) - L(1:k+1, 1:k-1)*(L(1:k+1, 1:k-1)'*L(1:k+1, k));
        L(1:k+1, k) = L(1:k+1, k)/norm(L(1:k+1, k));
        
        ql(1:k+1) = ql(1:k+1) - L(1:k+1, k)*(L(1:k+1, k)'*ql(1:k+1));
        ql(1:k+1) = ql(1:k+1) - L(1:k+1, k)*(L(1:k+1, k)'*ql(1:k+1));

        rr(k+1) = conj(ql(k+1)/ql(1));
      end
      ql = conj(ql)/conj(ql(1));
      %norm(rr-ql)
      f2(i,j) = rr.'*c;
    end
  end

  digits = 500;
  K = mp(K, digits);
  H = mp(H, digits);
  fmp = mp(z, digits);
  for i = 1:size(z,1),
    for j = 1:size(z,2),
      L = z(i,j)*K-H;  
      [Q, ~] = qr(L);
      fmp(i,j) = Q(:, end)'*c/conj(Q(1, end));
    end
  end

  fmp = double(fmp);
  
end
