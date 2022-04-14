## Author: Zoïs Moitier

## usage: [A, B] = _fct_AB_tp(nu, w, zeta, k_max)
##
## Get the value for the function A_nu(zeta) and B_nu(zeta) when zeta is close to 0
## using Taylor expansion see Table 2 in [Temme:1997]. The computation of the
## coefficients are done in the notebook `symbolic/comp_taylor_tp.ipynb`.
##

function [A, B] = _fct_AB_tp(nu, zeta, k_max)
  A_poly = {
    [1],
    [
       3.875233119897875110e-6
      -1.042960436782955530e-5
      -4.988652219516832030e-5
      -5.766301847639425080e-5
       1.540027672092350800e-4
       6.728876062209395540e-4
       7.064172724196895690e-4
      -1.463707463503144970e-3
      -4.444444444444444440e-3
    ],
    [
       5.240810645254741830e-5
      -1.044740083911794460e-4
      -3.513351434385566560e-4
      -2.698633097062688100e-4
       3.686607906143003560e-4
       6.937355413545889740e-4
    ],
    [
       2.341211902873772360e-4
      -2.478905546632297180e-4
      -3.542119714577438410e-4
    ]
  };

  B_poly = {
    [
      -1.068692482303864950e-7
      -5.790288733920437250e-7
      -7.726359892556073710e-7
       2.446810161235557900e-6
       1.301640251645853890e-5
       1.676987092017008960e-5
      -5.844357254566870890e-5
      -3.020604489992245090e-4
      -3.642848652199096040e-4
       1.625687162683573490e-3
       8.888888888888888890e-3
       1.799887214135533090e-2
    ],
    [
      -3.422607087563164680e-6
      -1.550546207672541230e-5
      -1.706623532653438110e-5
       4.105607390988507010e-5
       1.709853491354951200e-4
       1.690921480285995480e-4
      -3.820954145531625640e-4
      -1.394063079777365490e-3
      -1.492829532134291720e-3
    ],
    [
      -5.368400106135578370e-5
      -1.861483019310767530e-4
      -1.514935008908280400e-4
       2.528601609445752130e-4
       7.110486511670866890e-4
       5.522130767212927900e-4
    ],
    [
      -3.256754833263098310e-4
      -7.585627165879864240e-4
      -4.746177965599598080e-4
    ]
  };

  inv_nu2 = power(nu, -2);
  A = polyval(A_poly{k_max+1}, zeta);
  B = polyval(B_poly{k_max+1}, zeta);
  for k = k_max:(-1):1
    A = polyval(A_poly{k}, zeta) .+ inv_nu2 .* A;
    B = polyval(B_poly{k}, zeta) .+ inv_nu2 .* B;
  end
end