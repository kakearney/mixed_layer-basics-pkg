function B = nemparams2arrays(A)
%NEMPARAMS2ARRAYS Stores NEMURO parameters in arrays
%
% B = nemparams2arrays(A)
%
% This function rearranges the numerous parameters needed to run the NEMURO
% model (or its variants) into critter x 1 (for group-related) or critter x
% critter (for flux-related) critters.  This can make the parameters easier
% to access as a group, as opposed to the FORTRAN-style named parameters
% used in the original code.  In each array, the row and column dimensions
% correspond to the following NEMURO state variables:
%
% 1 = PS, small phytoplankton
% 2 = PL, large phytoplankton
% 3 = ZS, small zooplankton
% 4 = ZL, large zooplankton
% 5 = ZP, predatory zooplankton
% 6 = NO3, nitrate
% 7 = NH4, ammonium
% 8 = PON, particulate organic nitrogen
% 9 = DON, dissolved organic nitrogen
% 10 = SiOH4, silicic acid
% 11 = Opal, particulate silica
%
% Input variables:
%
%   A:  1 x 1 structure whose fields correspond to NEMURO parameter names,
%       and whose values are scalar parameter values.  See
%       nemuroParamSets.m for some examples.   
%
% Output variables:
%
%   B:  1 x 1 structure with the following fields:
%
%       alpha1:         Light Dissipation coefficient of sea water (/m)
% 
%       alpha2:         phytoplankton self-shading coefficient (l/molN/m) 
% 
%       usesteele:      logical scalar indicating whether to use Steele (1)
%                       or Platt (0) photosynthesis-irradiance curves
% 
%       RSiN:           Si/N ratio (molSi/molN)
% 
%       alpha:          nbsvx1 array, initial slope of
%                       photosynthesis-irradiance curve (Platt only)
%                       (/(ly/min)/s)   
% 
%       Iopt:           nbsvx1 array, optimal light intensity (Steele only)
%                       (ly/min)  
% 
%       Vmax:           nbsvx1 array, maximum photosynthetic rate @0degC
%                       (/s)  
% 
%       Kno3:           nbsvx1 array, half-saturation constant for uptake
%                       of NO3 (molN/l) 
% 
%       Knh4:           nbsvx1 array, hal-saturation constant for uptake of
%                       NH4 (molN/l) 
% 
%       Kgpp:           nbsvx1 array, temperature coefficient for
%                       photosynthetic rate (/degC) 
% 
%       pusai:          nbsvx1 array, ammonium inhibition coefficent
%                       (l/molN) 
% 
%       res0:           nbsvx1 array, respiration rate @0degC (/s)
% 
%       Kres:           nbsvx1 array, temperature coefficient for
%                       respiration (/degC) 
% 
%       mor0:           nbsvx1 array, mortality rate @0degC (l/(molN s))
% 
%       Kmor:           nbsvx1 array, temperature coefficient for mortality
%                       (/degC)  
% 
%       gamma:          nbsvx1 array, ratio of extracellular excretion to
%                       photosynthesis (no units) 
% 
%       lambda:         nbsvx1 array, Ivlev constant (l/molN)
% 
%       Kgra:           nbsvx1 array, temperature coefficient for grazing
%                       (/degC) 
% 
%       Ksi:            nbsvx1 array, half-saturation coefficient for
%                       silica (molSi/l) 
% 
%       alphaeg:        nbsvx1 array, assimilation effieciency (no units)
% 
%       beta:           nbsvx1 array, growth efficiency (no units)
% 
%       needsi:         nbsvx1 logical array, true if requires silica for
%                       growth 
% 
%       grmax:          nbsvxnbsv array, maximum rate of grazing by
%                       predator j on prey i @0degC (/s) 
% 
%       thresh:         nbsvxnbsv array, threshold value for grazing
%                       (molN/l) 
% 
%       grpusai:        nbsvxnbsv array, preference coefficient for prey i
%                       by predator j (l/molN) 
%
%       vdec:           nbsvxnbsv array, decomposition (or nitrification)
%                       rate of i into j @0degC (/s) 
%
%       Kdec:           nbsvxnbsv array, temperature coefficient for
%                       decomposition (/degC) 
%
%       inhibitedby:    nbsvxnbsv cell array, indices of prey groups
%                       preferred over prey i by predator j  
% 
%       settle:         nbsvx1 array, settling velocity (m/s)
%
%       egenosink:      nbsvx1 array, fraction of food egested by each
%                       grazer that does not sink (i.e. goes to DON instead
%                       of PON).  

% Number of state variables

nb = 11;

% Parameters that are scalars

B.alpha1 = A.alpha1;
B.alpha2 = A.alpha2;
B.usesteele = A.usesteele;
B.RSiN = A.RSiN;

% Parameters that are now nbsv x 1 arrays

[B.alpha, B.Iopt, B.Vmax, B.Kno3, B.Knh4, B.Kgpp, B.pusai, B.res0, ...
    B.Kres, B.mor0, B.Kmor, B.gamma, B.lambda, B.Kgra, B.Ksi, ...
    B.alphaeg, B.beta, B.needsi, B.settle, B.egenosink] = ...
    deal(zeros(nb,1));

B.alpha(1) = A.alphaS;
B.alpha(2) = A.alphaL;

B.Iopt(1) = A.IoptS;
B.Iopt(2) = A.IoptL;

B.Vmax(1) = A.VmaxS;
B.Vmax(2) = A.VmaxL;

B.Kno3(1) = A.KNO3S;
B.Kno3(2) = A.KNO3L;

B.Knh4(1) = A.KNH4S;
B.Knh4(2) = A.KNH4L;

B.Kgpp(1) = A.KGppS;
B.Kgpp(2) = A.KGppL;

B.pusai(1) = A.PusaiS;
B.pusai(2) = A.PusaiL;

B.res0(1) = A.ResPS0;
B.res0(2) = A.ResPL0;

B.Kres(1) = A.KResPS;
B.Kres(2) = A.KResPL;

B.mor0(1) = A.MorPS0;
B.mor0(2) = A.MorPL0;
B.mor0(3) = A.MorZS0;
B.mor0(4) = A.MorZL0;
B.mor0(5) = A.MorZP0;

B.Kmor(1) = A.KMorPS;
B.Kmor(2) = A.KMorPL;
B.Kmor(3) = A.KMorZS;
B.Kmor(4) = A.KMorZL;
B.Kmor(5) = A.KMorZP;

B.gamma(1) = A.GammaS;
B.gamma(2) = A.GammaL;

B.lambda(3) = A.LamS;
B.lambda(4) = A.LamL;
B.lambda(5) = A.LamP;

B.Kgra(3) = A.KGraS;
B.Kgra(4) = A.KGraL;
B.Kgra(5) = A.KGraP;

B.Ksi(2) = A.KSiL;
B.needsi = logical(B.needsi);
B.needsi(2) = true;

B.alphaeg(3) = A.AlphaZS;
B.alphaeg(4) = A.AlphaZL;
B.alphaeg(5) = A.AlphaZP;

B.beta(3) = A.BetaZS;
B.beta(4) = A.BetaZL;
B.beta(5) = A.BetaZP;

B.settle(nb-3) = A.setVPON;
B.settle(nb)   = A.setVOpal;

% Parameters that are now nbsv x nbsv (source x sink) arrays
% TODO: Check why I have GRmaxSps/pl for some param sets

[B.grmax, B.thresh, B.grpusai, B.vdec, B.Kdec] = deal(zeros(nb));
B.inhibitedby = cell(nb);

B.grmax(1,3) = A.GRmaxS;
B.grmax(2,4) = A.GRmaxLpl;
B.grmax(1,4) = A.GRmaxLps;
B.grmax(3,4) = A.GRmaxLzs;
B.grmax(2,5) = A.GRmaxPpl;
B.grmax(4,5) = A.GRmaxPzl;
B.grmax(3,5) = A.GRmaxPzs;

B.thresh(1,3) = A.PS2ZSstar;
B.thresh(2,4) = A.PL2ZLstar;
B.thresh(1,4) = A.PS2ZLstar;
B.thresh(3,4) = A.ZS2ZLstar;
B.thresh(2,5) = A.PL2ZPstar;
B.thresh(4,5) = A.ZL2ZPstar;
B.thresh(3,5) = A.ZS2ZPstar;

B.grpusai(2,5) = A.PusaiPL;
B.grpusai(3,5) = A.PusaiZS;

B.inhibitedby{2,5} = [3 4];
B.inhibitedby{3,5} = 4;

B.vdec(nb,  nb-1) = A.VO2S0;
B.vdec(nb-3,nb-2) = A.VP2D0;
B.vdec(nb-2,nb-4) = A.VD2N0;
B.vdec(nb-3,nb-4) = A.VP2N0;
B.vdec(nb-4,nb-5) = A.Nit0;    % nitrification uses same eq as decomp

B.Kdec(nb,  nb-1) = A.KO2S;
B.Kdec(nb-3,nb-2) = A.KP2D;
B.Kdec(nb-2,nb-4) = A.KD2N;
B.Kdec(nb-3,nb-4) = A.KP2N;
B.Kdec(nb-4,nb-5) = A.KNit;



