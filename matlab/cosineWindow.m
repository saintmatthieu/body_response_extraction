% cosineWindow Cosine window of order 0 to 3.
%     w_ = cosineWindow(N,P,cont) is the cosine window of size (N,1) and
%     order P. If cont is 0, the window is minimal-sidelobe, but with decay
%     of 6db/octave. If non-zero (default), the window is continuous at its
%     boundary to its (P-1)th derivative, which maximises the decay rate of
%     the sidelobes.
% 
%     [w_,a_] = cosineWindow(N,P,cont) returns the coefficient vector a_ of
%     the window, of length P+1 and with a_(1) such that sum(w_) equals N.

function [w_,a_] = cosineWindow(N,P,cont)

if nargin<4, cont = 1; end

switch P
case 0
    a_ = 1;
case 1
    if cont a_ = [1 1]';
    else a_ = [1 349/407]'; end
case 2
    if cont a_ = [21 25 4]'/21;
    else a_ = [1 1152/983 515/2792]'; end
case 3
    if cont a_ = [10 15 6 1]'/10;
    else a_ = [1 1667/1239 866/2305 161/5501]'; end
otherwise
    error('Window order P must be 0, 1, 2 or 3');
end

w_ = ones(N,1);
n_ = (0:N-1)';
for p=1:P
    w_(n_+1) = w_(n_+1)+(-1)^p*a_(p+1)*cos(p*2*pi*n_/N);
end

