function [x,success] = cplxsort(x)
%CPLXSORT   Sort vector by complex conjugate pairs.
%
% Used to sort poles in complex conjugate pairs if possible. Matlab's
% function cplxpair does the same but it will break with an error 
% if no pairing exists.

try
    x = cplxpair(x);
    success = 1;
catch
    success = 0;
end
end