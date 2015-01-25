function check = test_intf07()
  N   = 100;
  tol = 1e-12;

  A = gallery('grcar', N, 3);
  A = A/norm(full(A));
  B = gallery('ris', N);
  B = B/norm(full(B)) + 5*speye(N);
  b = ones(N, 1); 
  b(1)   = -10;
  b(end) = 5;

  xi = [inf, inf, -10i, 10i, -10+2i, -10-2i, ...
	inf, 200+45i, 200-45i, inf];
  
  inner_product = @(x, y) y'*(B*x);
  param.inner_product = inner_product;
  
  b = b/sqrt(inner_product(b,b));
  
  [V, K, H] = rat_krylov(A, B, b, xi, 'real', param);
  scale = norm(K)+norm(H);
  n1 = norm(A*V*K-B*V*H)/scale;
  n2 = norm(inner_product(V, V)-eye(size(V, 2)));
  n3 = ~isreal(K);

  [V, K, H] = rat_krylov(A, B, b, xi, param);
  scale = norm(K)+norm(H);
  n4 = norm(A*V*K-B*V*H)/scale;
  n5 = norm(inner_product(V, V)-eye(size(V, 2)));
  
  [Vr, Kr, Hr] = rat_krylov(B, A, b, H, K);
  n6 = norm(A*Vr*Hr-B*Vr*Kr)/scale;
  n7 = norm(inner_product(Vr, Vr)-eye(size(V, 2)))/scale;
  
  check = [n1 n2 n3 n4 n5 n6 n7] < tol;
end
