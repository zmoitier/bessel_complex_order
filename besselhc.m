function hnuz = besselhc(nu,l,z)
%BESSELHC Compute Hankel functions for order nu and argument z
% Both nu and z can be complex numbers.
% If nu is real uses the Matlab fonction besselh, if nu is non real use the
% uniform asymptotic expansion of Hankel functions in terms of Airy
% functions for large nu.
% (For details, see N.M. Temme, Numerical algorithms for uniform Airy-type
%  asymptotic expansions, Numerical Algorithms, vol. 15, p. 207-225, 1997,
%  section 2 and 3)
%
% Author : Z. Moitier, IRMAR, University of Rennes 1. April-June 2016.
% Last modification : 20 May 2016 (Zoïs Moitier)
% 
% usage :
%    jnuz = besselhc(nu,l,z)
%
% input parameters 
%    nu : [complex] order of the Bessel function.
%    l : [int] if k=1 return Hankel H1, if k=2 return Hankel H2.
%    z : [complex] argument of the Bessel function.
%
% output parameters
%    jnuz : [complex] value of Hankel function for order nu and argument z.
%    ierr : [int] If requested, ierr contains the following status
%                 information.
%

%
    if ((nargin<3)||(nargin>3))
        error('Invalid call to besselhc : wrong number of arguments, need 3 arguments');
    end
    %
    if ((length(nu)~=1)||(length(l)~=1)||(length(z)~=1))
        error('Invalid call to besselhc : nu, l and z must be scalar');
    end
    %
    if (and(l~=1,l~=2))
        error('Invalid call to besselhc : l=1 for Hankel H1 and l=2 for Hankel H2');
    end
    %
    if (imag(nu)==0)
        nu = real(nu);
        hnuz = besselh(nu,l,z);
    else
        icx = complex(0,1);
        jcx = exp(2*pi*icx/3);
        jbcx = conj(jcx);
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
            coeff1 = (4*dzeta/(1-y^2))^(0.25);
            if (l==1)
                % Compute H1_nu(z) by Formula (2.1) of [Temme:1997]
                coeffa = airy(0,jcx*nu^(2/3)*dzeta)/nu^(1/3);
                coeffb = jcx*airy(1,jcx*nu^(2/3)*dzeta)/nu^(5/3);
                hnuz = 2*exp(-icx*pi/3)*coeff1*(coeffa*suma+coeffb*sumb);
            else
                % Compute H2_nu(z) by Formula (2.1) of [Temme:1997]
                coeffa = airy(0,jbcx*nu^(2/3)*dzeta)/nu^(1/3);
                coeffb = jbcx*airy(1,jbcx*nu^(2/3)*dzeta)/nu^(5/3);
                hnuz = 2*exp(icx*pi/3)*coeff1*(coeffa*suma+coeffb*sumb);
            end
        else
            suma = 0;
            sumb = 0;
            [aktp,bktp] = abktp(dzeta);
            for k=0:3
                suma = suma + aktp(k+1)/nu^(2*k);
                sumb = sumb + bktp(k+1)/nu^(2*k);
            end
            if (l==1)
                % Compute H1_nu(z) by Formula (2.1) of [Temme:1997]
                coeffa = airy(0,jcx*nu^(2/3)*dzeta)/nu^(1/3);
                coeffb = jcx*airy(1,jcx*nu^(2/3)*dzeta)/nu^(5/3);
                hnuz = 2*exp(-icx*pi/3)*coeff1*(coeffa*suma+coeffb*sumb);
            else
                % Compute H2_nu(z) by Formula (2.1) of [Temme:1997]
                coeffa = airy(0,jbcx*nu^(2/3)*dzeta)/nu^(1/3);
                coeffb = jbcx*airy(1,jbcx*nu^(2/3)*dzeta)/nu^(5/3);
                hnuz = 2*exp(icx*pi/3)*coeff1*(coeffa*suma+coeffb*sumb);
            end
        end
    end
end
