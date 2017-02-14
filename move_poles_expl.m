function [KT, HT, QT, ZT] = move_poles_expl(K, H, xi)
% MOVE_POLES_EXPL    Changing the poles of the pencil (H, K).
%
% [KT, HT, QT, ZT] = move_poles_expl(K, H, xi) for (n+1)-by-n
% upper-Hessenberg matrices K and H and a vector xi of length k,
% with k < n+1, produces upper-Hessenberg matrices KT and HT and
% unitary matrices QT and ZT such that
%
%      KT = QT*K*ZT, HT = QT*H*ZT,
%
% and the first k poles of (KT, HT) are replaced by those specified
% by xi and pushed at the bottom.
%
% This algorithm is described in
%
% [1] M. Berljafa and S. G{\"u}ttel. Generalized rational Krylov
%     decompositions with an application to rational approximation,
%     SIAM J. Matrix Anal. Appl., 36(2):894--916, 2015.

  n = size(H, 2);
  k = length(xi);

  % Initial transformation to ensure upper-Hessenberg structure.
  if norm(abs(tril(H(3:end,:)))+abs(tril(K(3:end,:))), 'inf') ~= 0
    [H(2:end,:), K(2:end,:), QT, ZT] = qz(H(2:end,:), K(2:end,:));
    QT = blkdiag(1, QT);
    H(1, :) = H(1, :)*ZT;
    K(1, :) = K(1, :)*ZT;
  else
    QT = eye(n+1);
    ZT = eye(n);
  end

  %QT = eye(n+1);
  %ZT = eye(n);
  HT = H;
  KT = K;

  xi = xi(k:-1:1);

  for i = 1:k
    % Change the pole  HT(2, 1)/KT(2, 1).
    [s, c] = util_rat_givens(HT(1:2, 1), KT(1:2, 1), xi(i));
    KT(1:2, :) = [c -s;s' c]*KT(1:2, :);
    HT(1:2, :) = [c -s;s' c]*HT(1:2, :);
    QT(1:2, :) = [c -s;s' c]*QT(1:2, :);

    % Push the new HT(2, 1)/KT(2, 1) pole to the rear.
    index = ones(n-i+1, 1);
    index(1) = 0;

    [H, K, QS, ZS] = ordqz(HT(2:n-i+2, 1:n-i+1), ...
                           KT(2:n-i+2, 1:n-i+1), ...
                           eye(n-i+1), eye(n-i+1), ...
                           index);

    KT = [KT(1, 1:n-i+1)*ZS    KT(1,       n-i+2:end);
           K                QS*KT(2:n-i+2, n-i+2:end);
           KT(n-i+3:end, :)];
    HT = [HT(1, 1:n-i+1)*ZS    HT(1,       n-i+2:end);
           H                QS*HT(2:n-i+2, n-i+2:end);
           HT(n-i+3:end, :)];

    QT(2:n-i+2, :) = QS*QT(2:n-i+2, :);
    ZT(:, 1:n-i+1) = ZT(:, 1:n-i+1)*ZS;
  end
end