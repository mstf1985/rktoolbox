function [KT, HT, QT, ZT] = move_poles_impl(K, H, c) 
% MOVE_POLE_IMPL    Changing the poles of the pencil (H, K).
%
% [KT, HT, QT, ZT] = move_poles_impl(K, H, c) for (n+1)-by-n
% upper-Hessenberg matrices K and H and an (n+1)-by-1 vector c, 
% produces upper-Hessenberg matrices KT and HT and unitary
% matrices QT and ZT such that
% 
%      KT = QT*K*ZT, HT = QT*H*ZT,
%
% and the poles of (KT, HT) are replaced by those imposed by c.
%
% This algorithm is described in 
%
% [1] M. Berljafa and S. G\"{u}ttel. Generalized rational Krylov
%     decompositions with an application to rational approximation,
%     MIMS EPrint 2014.59, Manchester Institute for Mathematical
%     Sciences, The University of Manchester, UK, 2014. 

  c = c/norm(c);
  [u, beta] = house(c);
  HT = H - (beta*u)*(u'*H);
  KT = K - (beta*u)*(u'*K);
  
  [H, K, QT, ZT] = qz(HT(2:end, :), KT(2:end, :));
  QT = blkdiag(1, QT);
  HT = [HT(1, :)*ZT; H];
  KT = [KT(1, :)*ZT; K];
  
  QT = QT - (beta*(QT*u))*u';
end




function [v, beta] = house(x)  
%
%
  
  sigma(1) = x(2:end)'*x(2:end);
  v = [1; x(2:end)];
  if     sigma(1) == 0 && x(1) > 0, beta = 0;
  elseif sigma(1) == 0 && x(1) < 0, beta = -2;
  else
    mu = sqrt(sigma(1) + x(1)^2);
    if x(1) < 0, v(1) = x(1) - mu;
    else         v(1) = -sigma(1)/(x(1) + mu); end
    beta = 2*v(1)^2/(sigma(1) + v(1)^2);
    v = v/v(1);
  end
end
