function ynuz = besselyc(nu,z)
%BESSELYC Compute Bessel function of the second kind for order nu and
%         argument z.
% Both nu and z can be complex numbers.
% If nu is real uses the Matlab fonction bessely, if nu is non real use the
% uniform asymptotic expansion of Bessel functions in terms of Airy
% functions for large nu.
% (For details, see N.M. Temme, Numerical algorithms for uniform Airy-type
%  asymptotic expansions, Numerical Algorithms, vol. 15, p. 207-225, 1997,
%  section 2 and 3)
%
% Author : Z. Moitier, IRMAR, University of Rennes 1. April-June 2016.
% Last modification : 20 May 2016 (Zoïs Moitier)
% 
% usage :
%    jnuz = besselyc(nu,z)
%
% input parameters 
%    nu : [complex] order of the Bessel function.
%    z : [complex] argument of the Bessel function.
%
% output parameters
%    jnuz : [complex] value of Bessel function of the first second for
%                     order nu and argument z.
%

%
    if ((nargin<2)||(nargin>2))
        error('Invalid call to besseljc : wrong number of arguments, need 2 arguments');
    end
    %
    if ((length(z)~=1)||(length(nu)~=1))
        error('Invalid call to besseljc : nu and z must be scalar');
    end
    %
    if (imag(nu)==0)
        nu = real(nu);
        ynuz = bessely(nu,z);
    else
        y = z/nu;
        
        if (y==1)
            dzeta = 0;
            coeff1 = 2^(1/3);
            ak = 0;
            bk = 0;
            a2 = 1;
            b2 = 1;
        else
            transy = (1-y^2)^(0.5);
            % Formula (2.3) of [Temme:1997]
            if (abs(y) <= 1)
                dzeta = (1.5*(log((1+transy)/y)-transy))^(2/3);
            else
                dzeta = -(1.5*(sqrt(y^2-1.)-acos(1./y)))^(2/3);
            end
            coeff1 = (4.*dzeta/(1.-y^2))^(0.25);
            
            % Compute lambda et mu, cf Formula (2.6) of [Temme:1997]
            [lambda,mu] = lambdamu(5);
            
            % Compute Debye polynomial up to order 5, cf Formula (2.7) of [Temme:1997]
            polydebye = debye(5);
            suma = 0;
            sumb = 0;
            for k=0:1
                % Compute the coefficients given by Formula (2.5) of [Temme:1997]
                [ak,bk]=abk(k,dzeta,transy,lambda,mu,polydebye);
                suma = suma + ak/nu^(2*k);
                sumb = sumb + bk/nu^(2*k);
            end
            [a2,b2]=abk(2,dzeta,transy,lambda,mu,polydebye);
        end
        
        if (abs(a2) < abs(ak))&&(abs(b2) < abs(bk))
            suma = suma + a2/nu^4;
            sumb = sumb + b2/nu^4;
            % Compute Y_nu(z) by Formula (2.1) of [Temme:1997]
            coeff1 = (4*dzeta/(1-y^2))^(0.25);
            coeffa = airy(2,nu^(2/3)*dzeta)/nu^(1/3);
            coeffb = airy(3,nu^(2/3)*dzeta)/nu^(5/3);
            ynuz = -coeff1*(coeffa*suma+coeffb*sumb);
        else
            suma = 0;
            sumb = 0;
            [aktp,bktp] = abktp(dzeta);
            for k=0:3
                suma = suma + aktp(k+1)/nu^(2*k);
                sumb = sumb + bktp(k+1)/nu^(2*k);
            end
            % Compute Y_nu(z) by Formula (2.1) of [Temme:1997]
            coeffa = airy(2,nu^(2/3)*dzeta)/nu^(1/3);
            coeffb = airy(3,nu^(2/3)*dzeta)/nu^(5/3);
            ynuz = -coeff1*(coeffa*suma+coeffb*sumb);
        end
    end
end
