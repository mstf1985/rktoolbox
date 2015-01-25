function check = test_rkfit()

A = [ 3 , 1 ; 1 2-3i ];
F = expm(A);
v = [ -1i;1 ];
poles = 1;
[poles,ratfun,misfit] = rkfit(F,A,v,poles);
check(1) = (misfit/norm(F*v)) < 1e-14;
check(2) = abs(ratfun(poles+10*eps)) > 1e14;
