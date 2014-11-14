function [b1, sh, sm, Kq, Km, Kh] = mycoef(nz, sh, gh, sm, len, q2, kmol)
%MYCOEF Calculate mixing coefficients
%
% [b1, sh, sm, Kq, Km, Kh] = mycoef(nz, sh, gh, sm, len, q2)
%
% This routine calculates the mixing coefficients for turbulence, scalars
% and momentum.  A typical value for molecular diffusion (kmol) is added to
% the calculated values to serve as a lowe limit for diffusion.

% Charlie Stock
% Updated by Kelly Kearney, Aug 2008



% Constants relating various length-scales to master length scale
a1 = 0.92;
b1 = 16.6;
a2 = 0.74;
b2 = 10.1;
c1 = 0.08;

% solve for new mixing coefficients
coef1 = a2*(1 - 6*a1/b1);
coef2 = 3*a2*b2 + 18*a1*a2;
coef3 = a1*(1-3*c1-6*a1/b1);
coef4 = 18*a1*a1+9*a1*a2;
coef5 = 9*a1*a2;

% solve for new coefficients, if no wave-breaking, the first index 
% should be 2.  If wave-breaking, the first index should be 1 

sh(1:nz) = coef1./(1-coef2*gh(1:nz));
sm(1:nz) = (coef3 + sh(1:nz).*coef4.*gh(1:nz))./(1-coef5*gh(1:nz));
Kn(1:nz) = len(1:nz).*sqrt(q2(1:nz));

[Kq, Km, Kh] = deal(zeros(size(sh)));

Kq(1:nz) = 0.41*sh(1:nz).*len(1:nz).*sqrt(q2(1:nz)) + kmol;
Km(1:nz) = sm(1:nz).*len(1:nz).*sqrt(q2(1:nz)) + kmol;
Kh(1:nz) = sh(1:nz).*len(1:nz).*sqrt(q2(1:nz)) + kmol;
