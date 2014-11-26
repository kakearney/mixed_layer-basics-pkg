function ps = photosynthesis(nutrients, phyto, kn, z, irr, temp, dz)
%PHOTOSYNTHESIS Calculate growth due to photosynthesis
%
% ps = photosynthesis(nutrients, phyto, kn, z, irr, temp, dz)
%
% Using photosynthesis model from Sarmiento and Gruber 2006, Chapter 4.
% Net uptake and assimilation of nitrogen is defined as Vp(T)*gammap(I,N),
% where Vp(T) is the maximum temperature-dependant growth rate (time^-1),
% and gammap(I,N) is a value between 0 and 1 that describes the limitation
% of growth due to light and nutrients.
% 
% Input variables:
%
%   phyto:      ndepth x np, concentration of producers (mmol N m^-3)
%
%   nutrients:  ndepth x 1, concentration of nutrients (mmol N m^-3)
%
%   kn:         1 x np, half-saturation constants for each producer (mmol N
%               m^-3)
%
%   z:          ndepth x 1, depths (m, positive down)
%
%   irr:        mean irradiance at surface (W m^-2)
%
%   dz:         ndepth x 1 or 1 x 1 array, thickness of each layer (m)
%
% Output variables:
%
%   ps:         ndepth x np, nutrient uptake rate (or rate of accumulation
%               of biomass due to photosynthesis) per unit biomass of
%               producer (s^-1)


% Temperature dependance of maximum photosynthetic rate according to Eppley
% 1972 (per depth)

a = 0.6;        % d^-1
a = a ./ 86400; % s^-1
b = 1.066;      % no units
c = 1;          % (deg C)^-1

vmax = a.*b.^(c.*temp); % s^-1

% Influence of nutrients on growth (per depth layer and producer group)

gamman = bsxfun(@rdivide, nutrients, bsxfun(@plus, kn, nutrients));

% Influence of light on growth, using Platt and Jasby 1976

kw = 0.04;  % m^-1, clear water attenuation coefficient
kp = 0.03;  % m^-1 (mmol m^-3)^-1, phytoplankton effect on light
phyto2 = sum(phyto, 2); % All producers treated same for light absorption
kpp = cumsum(kp .* phyto2 .* dz)./z;

irrz = irr .* exp(-(kw + kpp).*z);  % depth-dependant irradiance

alpha = 0.025;          % d^-1(W m^-2)^-1
alpha = alpha ./ 86400; % s^-1(W m^-2)^-1 

ik = vmax ./ alpha;

gammai = irrz ./ sqrt(ik.^2 + irrz.^2);

% Light and nutrient limitation

gammap = bsxfun(@times, gamman, gammai);

% Uptake of nutrient due to photosynthesis

ps = bsxfun(@times, vmax, gammap);
