function check = test_intf06()
  N   = 100;
  tol = 1e-10;

  A = gallery('riemann', N);
  A = A/norm(full(A));
  B = gallery('ris', N);
  B = B/norm(full(B)) + 5*speye(N);
  b = ones(N, 1); 
  b(1)   = 10;
  b(end) = 5;

  xi = zeros(1, 10);
  xi(1) = 1+1i; xi(2) = 1-1i;
  inner_product = @(x, y) y'*(B*x);
  param.inner_product = inner_product;
  
  b = b/sqrt(inner_product(b,b));
  
  [V, K, H] = rat_krylov(A, b, xi, param);
  scale = norm(K)+norm(H);
  n1 = norm(A*V*K-V*H)/scale;
  n2 = norm(inner_product(V, V)-eye(size(V, 2)));

  [Vr, Kr, Hr] = rat_krylov(A, b, K, H);
  n3 = norm(A*Vr*K-Vr*H)/scale;
  n4 = norm(V-Vr)/scale/scale;
  
  [Vr, Kr, Hr] = rat_krylov(B, A\b, K, H);
  n5 = norm(B*Vr*K-Vr*H)/scale;
   
  [Vr, Kr, Hr] = rat_krylov(B*A, B, b, K, H);
  n6 = norm(A*Vr*Kr-Vr*Hr)/scale;
  n7 = norm(Vr-V)/scale/scale;
  
  check = [n1 n2 n3 n4 n5 n6 n7] < tol;
end
