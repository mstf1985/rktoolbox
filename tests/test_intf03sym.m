function check = test_intf03sym()
  N   = 200;
  tol = 1e-10;

  A = speye(N);
  B = gallery('tridiag', N) + speye(N);
  b = ones(N, 1); 
  
  param.linear_solver = @(M, x) (M + randn(size(M))*1e-14)\x;
  param.inner_product = @(x, y) y'*(B*x);
  
  b = b/sqrt(param.inner_product(b, b));
  
  xi = repmat([-1, -5, -10, -20, -25], 1, 8);

  [V, K, H] = rat_krylov(A, B, b, xi, 'real', param);
  scale = norm(K)+norm(H);
  n1 = norm(A*V*K-B*V*H)/scale;
  n2 = norm(param.inner_product(V, V)-eye(size(V, 2)));
  
  [Vr, Kr, Hr] = rat_krylov(A, B, b, K, H, param);
  n3 = norm(A*Vr*Kr-B*Vr*Hr)/scale;
  n4 = norm(V*param.inner_product(Vr, V)-Vr)/scale;
  
  C = gallery('toeppen', N);
  scale = scale*norm(full(B));
  [Vr, Kr, Hr] = rat_krylov(C*A, C*B, b, K, H, 'real');
  n5 = norm(C*A*Vr*K-C*B*Vr*H)/scale;
  n6 = norm(V-Vr)/scale;
  
  check = [n1 n2 n3 n4 n5 n6] < tol;
end
