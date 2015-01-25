function check = test_intf05expl()
  N   = 100;
  tol = 1e-12;

  A = gallery('tridiag', N);
  B = gallery('dramadah', N, 3);
  b = ones(N, 1);
  b = b/norm(b);
  
  xi = [123, 17+2i, 17-2i, 22+10i, 22-10i];
  
  [V, K, H] = rat_krylov(A, B, b, xi, 'real');
  scale = norm(K)+norm(H);
  n1 = norm(A*V*K-B*V*H)/scale;
  n2 = norm(V'*V-eye(size(V, 2)));

  [KT, HT, QT, ZT] = move_poles_expl(K, H, inf*xi);
  VT = V*QT';
  n3 = norm(A*VT*KT-B*VT*HT)/scale;
      
  [Vr, Kr, Hr] = rat_krylov(A, B, VT(:, 1), KT, HT);
  n4 = norm(A*Vr*KT-B*Vr*HT)/scale;
  n5 = norm(V*(V'*Vr) - Vr);
  check = [n1 n2 n3 n4 n5] < tol;
end
