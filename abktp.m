function [aktp,bktp] = abktp(dzeta)
%ABKTP Compute the k-ieme coefficient of the a_k and b_k sequense for dzeta
%      close to the turnig point
% For details, see N.M. Temme, Numerical algorithms for uniform Airy-type
% asymptotic expansions, Numerical Algorithms, vol. 15, p. 207-225, 1997,
% section 3
%
% Author : Z. Moitier, IRMAR, University of Rennes 1. April-June 2016.
% Last modification : 24 May 2016 (Zoïs Moitier)
% 
% usage :
%  [aktp,bktp] = abktp(dzeta)
%
% input parameters 
%    dzeta : [complex]
%
% output parameters
%    ak : [complex] coefficient a_k
%    bk : [complex] coefficient b_k
%

%
    eta = dzeta*0.79370052598409974;
    
    poly = zeros(1,6);
    aktp = zeros(1,4);
    bktp = zeros(1,4);
    
    aktp(1) = 1;
   	
    poly(2) = 3.8806265629795042e-4;
    poly(3) = 1.3457752124418791e-3;
    poly(4) = 1.1213675213675214e-3;
    poly(5) = -1.8441558441558442e-3;
    poly(6) = -4.4444444444444444e-3;
    aktp(2) = polyval(poly,eta);
    
    poly(2) = 0;
    poly(3) = -7.0267028687711331e-4;
    poly(4) = -4.2838130171535112e-4;
    poly(5) = 4.6448349036584331e-4;
    poly(6) = 6.9373554135458897e-4;
    aktp(3) = polyval(poly,eta);
    
    poly(3) = 0;
    poly(4) = 3.7164422375022963e-4;
    poly(5) = -3.1232252789031883e-4;
    poly(6) = -3.5421197145774384e-4;
    aktp(4) = polyval(poly,eta);
    
    poly(1) = -1.8554677707954858e-4;
    poly(2) = -7.6114463606963947e-4;
    poly(3) = -7.2856973043981921e-4;
    poly(4) = 2.580617512215102e-3;
    poly(5) = 1.1199298221287761e-2;
    poly(6) = 1.7998872141355331e-2;
    bktp(1) = polyval(poly,eta);
    
    poly(1) = 0;
    poly(2) = 4.3085608119886891e-4;
    poly(3) = 3.381842960571991e-4;
    poly(4) = -6.0653866301391553e-4;
    poly(5) = -1.7564094190927787e-3;
    poly(6) = -1.4928295321342917e-3;
    bktp(2) = polyval(poly,eta);
    
    poly(2) = 0;
    poly(3) = -3.0298700178165608e-4;
    poly(4) = 4.0139048548426692e-4;
    poly(5) = 8.9586516310476929e-4;
    poly(6) = 5.5221307672129279e-4;
    bktp(3) = polyval(poly,eta);
    
    poly(3) = 0;
    poly(4) = -5.1697760483243603e-4;
    poly(5) = -9.5572913429464297e-4;
    poly(6) = -4.7461779655995981e-4;
    bktp(4) = polyval(poly,eta);
end
