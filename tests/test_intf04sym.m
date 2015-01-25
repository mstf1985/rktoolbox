function check = test_intf04sym()
  N   = 400;
  tol = 1e-10;

  A = gallery('tridiag', N, 1, 100, 0.45);
  B = speye(N);
  b = ones(N, 1); 
  
  b = b/norm(b);
  
  xi = rand(1, 5) + 2*rand(1, 5)*1i;
  xi = xi - 100;
  xi = [xi conj(xi)];
  xi = cplxsort(xi);
  xi = repmat(xi, 1, 3);
  xi = [xi -50 inf];
  
  [V, K, H] = rat_krylov(A, B, b, xi, 'real');
  scale = norm(K)+norm(H);
  n1 = norm(A*V*K-B*V*H)/scale;
  n2 = norm(V'*V-eye(size(V, 2)));
  
  [Vr, Kr, Hr] = rat_krylov(A, B, b, K, H);
  n3 = norm(A*Vr*K-B*Vr*H)/scale;
  n4 = norm(V*(V'*Vr)-Vr)/scale;
 
  [Vr, Kr, Hr] = rat_krylov(B\A, b, K, H);
  n5 = norm(A*Vr*K-B*Vr*H)/scale;
  n6 = norm(V-Vr)/scale;
    
  [Vp, Kp, Hp] = rat_krylov(A, V(:, end), inf*xi);
  n7 = norm(A*Vp*Kp-B*Vp*Hp)/scale;
  n8 = norm(V*(V'*Vp)-Vp);
    
  check = [n1 n2 n3 n4 n5 n6 n7] < tol;
end
