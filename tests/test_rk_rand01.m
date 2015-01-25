function check = test_rk_rand01(varargin)
  n = 19;
  N = n*n;
  
  tol = 1e-13*N;
  
  b = ones(N, 1)/n; 
  b(2) = 2;
    
  nA  = 5;
  nxi = 12;
  
  check = zeros(nA, nxi);
  
  A  = cell(1, nA);
  B  = cell(1, nA);
  xi = cell(1, nxi);
  
  B{2} = speye(N);
  B{3} = 2*speye(N);
  B{1} = gallery('tridiag', N)-100*speye(N);
  B{4} = speye(N);
  B{5} = 1i*speye(N);
  
  A{1} = gallery('tridiag', N);
  A{2} = kron(gallery('tridiag', sqrt(N)), ...
              gallery('tridiag', sqrt(N)));
  A{3} = gallery('grcar', N, 3);
  A{4} = gallery('toeppen', N);
  A{5} = gallery('toeppen', N) + 1i*gallery('tridiag', N);
  
  xi{1}  = [-10, 10, inf];
  xi{2}  = [10, inf, 1i, -1i];
  xi{3}  = 2*[1+2i, 1-2i, 1+3i, 1-3i, 1-2i, 1+2i];
  xi{4}  = [-3, -4, -6, -12, -4, -4];
  xi{5}  = [4-2i, 4+2i];
  xi{6}  = [4-2i, 4+2i, 4+2i, 4-2i, 4-2i, 4+2i, 4+2i, 4-2i];
  xi{7}  = [-3, -4, -6, -12, -4, -4, 4-2i, 4+2i, 4-2i, 4+2i, 4+2i, 4-2i];
  xi{8}  = [0, 100-100i, 100+100i, 100+100i, 100-100i, 100-100i, ...
            100+100i, 0, inf, 0];
  xi{9}  = [-2, -3, -20];
  xi{10} = [-1i, 1i, -10i, 10i];
  xi{11} = [inf, inf, inf, inf, inf];
  xi{12} = [0, 0, 0, 0, 0, 0, 0, 0]-5;
  
  for i = 1:nA
    for j = 1:nxi
      if isreal(A{i}) && isreal(B{i})
        [V, K, H] = rat_krylov(A{i}, B{i}, b, xi{j}, 'real');   
        
        n1 = norm(A{i}*V*K-B{i}*V*H)/(norm(K)+norm(H));
        n2 = norm(V'*V-eye(size(V, 2)));     
        
        check(i, j) = n1 < tol && n2 < tol;
      else
        check(i, j) = true;
      end
      
      [V, K, H] = rat_krylov(A{i}, B{i}, b, xi{j});
      
      n1 = norm(A{i}*V*K-B{i}*V*H)/(norm(K)+norm(H));
      n2 = norm(V'*V-eye(size(V, 2)));
      
      check(i, j) = check(i, j) && n1 < tol && n2 < tol;
    end
  end

end