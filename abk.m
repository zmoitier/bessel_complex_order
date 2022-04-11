function [ak,bk] = abk(k,dzeta,transy,lambda,mu,polydebye)
%ABK Compute the k-ieme coefficient of the a_k and b_k sequense
% For details, see N.M. Temme, Numerical algorithms for uniform Airy-type
% asymptotic expansions, Numerical Algorithms, vol. 15, p. 207-225, 1997,
% section 2
%
% Author : Z. Moitier, IRMAR, University of Rennes 1. April-June 2016.
% Last modification : 17 May 2016 (Zoïs Moitier)
% 
% usage :
%  [ak,bk] = abk(k,dzeta,transy,lambda,mu,polydebye)
%
% input parameters 
%    k : [int] k-ieme coefficient
%    dzeta : [complex]
%    transy : [complex] non zeros
%    lambda : [(NbPts x 1) array] lambda sequence
%    mu : [(NbPts x 1) array] mu sequence
%    polydebye : [cell-array] the list of coefficients of Debye polynoms
%
% output parameters
%    ak : [complex] coefficient a_k
%    bk : [complex] coefficient b_k
%

%
    z = 1/transy;
    % initialization
    ak = mu(1)*polyval(polydebye{2*k+1},z);
    bk = lambda(1)*polyval(polydebye{2*k+2},z);
    % boucle to compute the sum
    for s=1:(2*k)
        transd = dzeta^(-1.5*s);
        ak = ak + mu(s+1)*transd*polyval(polydebye{2*k-(s-1)},z);
        bk = bk + lambda(s+1)*transd*polyval(polydebye{2*k-(s-1)+1},z);
    end
    bk = -dzeta^(-0.5)*(bk+lambda(2*k+2)*dzeta^(-1.5*(2*k+1))*polyval(polydebye{1},z));
end
