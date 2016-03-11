function Qo = clearsky(start_date,t,Lat)
%CLEARSKY Calculate clear-sky irradiance
%
% Qo = clearsky(start_date,t,Lat)
%
% This function calculates the clear-sky irradiance using the formula of
% Seckel and Beaudry (1973) as reported in Reed (1977).  This is used in
% the calculation of net longwave heat flux.
%
% Input variables:
%
%   start_date: 1 x 6 date vector of time corresponding to t = 0
%
%   t:          vector of time values to calculate clear-sky irradiance
%               for, in seconds since start_date 
%
%   Lat:        scalar, latitude of location of interest
%
% Output variables:
%
%   Qo:         vector same length as t, clear-sky irradiance (W m^-2)



tyear = datenum(start_date) + (t/86400) - ...
        datenum([start_date(1) 1 1 0 0 0]);
phi = pi*((tyear-21)*(360/365))/180;

if ((Lat >= 40) & (Lat <= 60))
   A0 = 342.61 - 1.97*Lat - 0.018*(Lat^2);
   A1 = 52.08 - 5.86*Lat + 0.043*(Lat^2);
   B1 = -4.80 + 2.46*Lat - 0.017*(Lat^2);
   A2 = 1.08 - 0.47*Lat + 0.011*(Lat^2);
   B2 = -38.79 + 2.43*Lat - 0.034*(Lat^2);
elseif ((Lat >= -20) & (Lat < 40))
   A0 = -15.82 + 326.87*cos(pi*Lat/180);
   A1 = 9.63 + 192.44*cos(pi*(Lat+90)/180);
   B1 = -3.27 + 108.70*sin(pi*Lat/180);
   A2 = -0.64 + 7.80*(sin(pi*(Lat-45)/180))^2;
   B2 = -0.50 + 14.42*(cos(pi*(Lat-5)/180))^2;
else
   disp(sprintf('The Seckel and Beaudry formula for clear sky\n'));
   disp(sprintf('is only valid for latitudes between 20S and 60N.\n'));
   disp(sprintf('At other latitudes, the cloud correction factors\n'));
   disp(sprintf('for the calculation of Qlw will be incorrect.'));
end

Qo = A0+A1*cos(phi)+B1*sin(phi)+A2*cos(2*phi)+B2*sin(2*phi);