function check = test_realopt2(varargin)

  n = 14;
  N = n*n;

  tol = 1e-10;

  b = ones(N, 1)/n;
  b = randn(N, 1); b = b/norm(b);

  nA  = 6; nAr = 3;
  nxi = 10;

  A  = cell(1, nA);
  Ar = cell(1, nAr);
  xi = cell(1, nxi);

  check = zeros(nA, nxi, nAr);

  A{1} = gallery('tridiag', N);
  A{2} = kron(gallery('tridiag', sqrt(N)), ...
              gallery('tridiag', sqrt(N)));
  A{3} = gallery('grcar', N, 3);
  A{4} = gallery('grcar', N, 4);
  A{5} = gallery('toeppen', N);
  A{6} = gallery('toeppen', N) + 1i*gallery('tridiag', N);;

  Ar{1} = gallery('tridiag', N) + speye(N);
  Ar{2} = gallery('tridiag', N) + 10i*speye(N);
  Ar{3} = (16-1i)*kron(gallery('tridiag', sqrt(N)), ...
                       gallery('tridiag', sqrt(N)));

  xi{1}  = [-10, 10, inf];
  xi{2}  = [2+1i, 2-1i];
  xi{3}  = 2*[10, 1+3i, 1-3i, 1-2i, 1+2i];
  xi{4}  = [-3, -4, -6, -12, -4, -4];
  xi{5}  = [4-2i, 4+2i];
  xi{6}  = [4-2i, 4+2i, 4+2i, 4-2i, 4-2i, 4+2i, 4+2i, 4-2i];
  xi{7}  = [-3, -4, -6, 4-2i, 4+2i, 4-2i, 4+2i];
  xi{8}  = [0, 100-100i, 100+100i];
  xi{9}  = [-2, -3, -20, 0];
  xi{10} = [-1i, 1i, -1i, 1i, inf];

  for i = 1:nA
    for j = 1:nxi
      if isreal(A)
	[V, K, H]    = rat_krylov(A{i}, b, xi{j}, 'real');
	n1 = norm(A{i}*V*K-V*H)/(norm(K)+norm(H));
	n2 = norm(V'*V-eye(size(V, 2)));
	for k = 1:nAr
	  [Vr, Kr, Hr] = rat_krylov(Ar{k}, b/norm(b), K, H);
	  n3 = norm(Ar{k}*Vr*Kr-Vr*Hr)/(norm(Kr)+norm(Hr));
	  check(i, j, k) = n1 < tol && n2 < tol && n3 < tol;
	end
      else
	check(i, j, :) = ones(nAr, 1);
      end
      [V, K, H]    = rat_krylov(A{i}, b, xi{j});
      n1 = norm(A{i}*V*K-V*H)/(norm(K)+norm(H));
      n2 = norm(V'*V-eye(size(V, 2)));
      for k = 1:nAr
        [Vr, Kr, Hr] = rat_krylov(Ar{k}, b/norm(b), K, H);
        n3 = norm(Ar{k}*Vr*Kr-Vr*Hr)/(norm(Kr)+norm(Hr));
        check(i, j, k) = check(i, j, k) && n1 < tol && n2 < tol && ...
            n3 < tol;
      end
    end
  end

end