function [KT, HT, QT, ZT] = move_poles_impl(K, H, c, flag)
% MOVE_POLES_IMPL    Changing the poles of the pencil (H, K).
%
% [KT, HT, QT, ZT] = move_poles_impl(K, H, c, flag) for
% (n+1)-by-n (quasi-)upper-Hessenberg matrices K and H and an
% (n+1)-by-1 vector c, produces (quasi-)upper-Hessenberg matrices
% KT and HT and unitary matrices QT and ZT such that  
%
%      KT = QT*K*ZT, HT = QT*H*ZT,
%
% and the poles of (KT, HT) are replaced by those imposed by c.
%
% The parameter flag is optional, and can be passed as flag =
% 'real' to obtain quasi-upper-Hessenberg matrices KT and HT,
% provided that K, H, and c are real-valued. 
% 
% This algorithm is described in
%
% [1] M. Berljafa and S. G{\"u}ttel. Generalized rational Krylov
%     decompositions with an application to rational approximation,
%     SIAM J. Matrix Anal. Appl., 36(2):894--916, 2015.

  if nargin == 3, flag = 'complex'; end

  [u, beta, gamma] = util_householder(c);
  
  HT = H - (conj(beta)*u)*(u'*H);
  KT = K - (conj(beta)*u)*(u'*K);
   
  [KT, HT, QT, ZT] = util_recover_rad(KT, HT, flag);

  QT = QT - (conj(beta)*(QT*u))*u';
  
  HT = gamma*HT;
  KT = gamma*KT;
  QT = gamma*QT;    
end
