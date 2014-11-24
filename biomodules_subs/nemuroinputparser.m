function [A,B] = nemuroinputparser(varargin)
%NEMUROINPUTPARSER Checks and returns user input and defaults for NEMURO
%
% Params = nemuroinputparser('param1', val1, 'param2', val2, ...)
% Params = nemuroinputparser(setname)
% Params = nemuroinputparser(setname, 'param1', val1, 'param2', val2, ...)
% [Params, Extra] = nemuroinputparser(...)
% 
% Input variables:
%
%   setname:    string corresponding to one of the set names in
%               nemuroParamSets.mat.  These indicate various sets of NEMURO
%               parameters that have been published in the literature.  Any
%               unspecified variables (see below) will be set to a default
%               value from the specified set.  If no set is specified, the
%               'NEMURO Version 1.f90' set, from the original source code
%               representing the A7 location, will be used.
%                 'Eslinger et al. Simulation Parameter'
%                 'Eslinger et al. Station P'
%                 'Eslinger et al. A7'
%                 'Eslinger et al. Bering'
%                 'NEMURO Version 1.f90'
%                 'Fujii et al. Papa'
%                 'Fujii et al. A7'
%                 'Fujii et al. KNOT'
%                 'Kishi et al. A7'
%
% Optional input variables:
%
%     All variables should be scalars.  Values can be entered as
%     parameter/value pairs or as fields in a structure.  Default values
%     come from parameter set indicated by setname.
%
%     alpha1:       Light Dissipation coefficient of sea water  [/m]
%     alpha2:       PS+PL Selfshading coefficientS+PL           [l/molN/m]
%     IoptS:        PS Optimum Light Intensity  S               [ly/min]
%     IoptL:        PL Optimum Light Intensity                  [ly/min]
%     LLN:          Number of sublayer for calculating of Lfc
%     VmaxS:        PS Maximum Photosynthetic rate @0degC       [/s]
%     KNO3S:        PS Half satuation constant for Nitrate      [molN/l]
%     KNH4S:        PS Half satuation constant for Ammonium     [molN/l]
%     PusaiS:       PS Ammonium Inhibition Coefficient          [l/molN]
%     KGppS:        PS Temp. Coeff. for Photosynthetic Rate     [/degC]
%     MorPS0:       PS Mortality Rate @0degC                    [/s]
%     KMorPS:       PS Temp. Coeff. for Mortality               [/degC]
%     ResPS0:       PS Respiration Rate at @0degC               [/s]
%     KResPS:       PS Temp. Coeff. for Respiration             [/degC]
%     GammaS:       PS Ratio of Extracell. Excret. to Photo.    [(nodim)]
%     VmaxL:        PL Maximum Photosynthetic rate @0degC       [/s]
%     KNO3L:        PL Half satuation constant for Nitrate      [molN/l]
%     KNH4L:        PL Half satuation constant for Ammonium     [molN/l]
%     KSiL:         PL Half satuation constant for Silicate     [molSi/l]
%     PusaiL:       PL Ammonium Inhibition Coefficient          [l/molN]
%     KGppL:        PL Temp. Coeff. for Photosynthetic Rate     [/degC]
%     MorPL0:       PL Mortality Rate @0degC                    [/s]
%     KMorPL:       PL Temp. Coeff. for Mortality               [/degC]
%     ResPL0:       PL Respiration Rate at @0degC               [/s]
%     KResPL:       PL Temp. Coeff. for Respiration             [/degC]
%     GammaL:       PL Ratio of Extracell. Excret. to Photo.    [(nodim)]
%     GRmaxS:       ZS Maximum Rate of Grazing PS @0degC        [/s]
%     KGraS:        ZS Temp. Coeff. for Grazing                 [/degC]
%     LamS:         ZS Ivlev constant                           [l/molN]
%     PS2ZSstar:    ZS Threshold Value for Grazing PS           [molN/l]
%     AlphaZS:      ZS Assimilation Efficiency                  [(nodim)]
%     BetaZS:       ZS Growth Efficiency                        [(nodim)]
%     MorZS0:       ZS Mortality Rate @0degC                    [/s]
%     KMorZS:       ZS Temp. Coeff. for Mortality               [/degC]
%     GRmaxLps:     ZL Maximum Rate of Grazing PS @0degC        [/s]
%     GRmaxLpl:     ZL Maximum Rate of Grazing PL @0degC        [/s]
%     GRmaxLzs:     ZL Maximum Rate of Grazing ZS @0degC        [/s]
%     KGraL:        ZL Temp. Coeff. for Grazing                 [/degC]
%     LamL:         ZL Ivlev constant                           [l/molN]
%     PS2ZLstar:    ZL Threshold Value for Grazing PS           [molN/l]
%     PL2ZLstar:    ZL Threshold Value for Grazing PL           [molN/l]
%     ZS2ZLstar:    ZL Threshold Value for Grazing ZS           [molN/l]
%     AlphaZL:      ZL Assimilation Efficiency                  [(nodim)]
%     BetaZL:       ZL Growth Efficiency                        [(nodim)]
%     MorZL0:       ZL Mortality Rate @0degC                    [/s]
%     KMorZL:       ZL Temp. Coeff. for Mortality               [/degC]
%     GRmaxPpl:     ZP Maximum rate of grazing PL @0degC        [/s]
%     GRmaxPzs:     ZP Maximum rate of grazing ZS @0degC        [/s]
%     GRmaxPzl:     ZP Maximum rate of grazing ZL @0degC        [/s]
%     KGraP:        ZP Temp. Coeff. for grazing                 [/degC]
%     LamP:         ZP Ivlev constant                           [l/molN]
%     PL2ZPstar:    ZP Threshold Value for Grazing PL           [molN/l]
%     ZS2ZPstar:    ZP Threshold Value for Grazing ZS           [molN/l]
%     ZL2ZPstar:    ZP Threshold Value for Grazing ZL           [molN/l]
%     PusaiPL:      ZP Preference Coeff. for PL                 [l/molN]
%     PusaiZS:      ZP Preference Coeff. for ZS                 [l/molN]
%     AlphaZP:      ZP Assimilation Efficiency                  [(nodim)]
%     BetaZP:       ZP Growth Efficiency                        [(nodim)]
%     MorZP0:       ZP Mortality Rate @0degC                    [/s]
%     KMorZP:       ZP Temp. Coeff. for Mortality               [/degC]
%     Nit0:         NH4 Nitrification Rate @0degC               [/s]
%     KNit:         NH4 Temp. coefficient for Nitrification     [/degC]
%     VP2N0:        PON eecomp. Rate to Ammonium @0degC         [/s]
%     KP2N:         PON Temp. Coeff. for Decomp. to Ammon.      [/degC]
%     VP2D0:        PON Decomp. Rate to DON @0degC              [/s]
%     KP2D:         PON Temp. Coeff. for Decomp. to DON         [/degC]
%     VD2N0:        DON Decomp. Rate to Ammonium @0degC         [/s]
%     KD2N:         DON Temp. Coeff. for Decomp. to Ammon.      [/degC]
%     VO2S0:        Opal Decomp. Rate to Silicate @0degC        [/s]
%     KO2S:         Opal Temp. Coeff. for Decomp.to Silicate    [/degC]
%     RSiN:         Si/N ratio                                  [molSi/molN]
%     RCN:          C/N ratio                                   [molC/molN]
%     setVPON:      Settling velocity of PON                    [m/s]
%     setVOpal:     Settling velocity of Opal                   [m/s]
%     TNO3d:        Nitrate Concentraion in the Deep Layer      [molN/l]
%     TSiOH4d:      Silicate Concentraion in the Deep Layer     [molSi/l]
%     usesteele:    logical scalar indicating whether to use Steele (1) or
%                   Platt (0) light curves
%
% Output variables:
%
%     Params:       1 x 1 structure with fields corresponding to variables
%                   listed above, with default values if user did not
%                   provide a value via input.
%
%     Extra:        1 x 1 structure including any parameter-value pairs
%                   included by the user that are not included in the above
%                   list (used by nemuroflexinput)

% Copyright 2009 Kelly Kearney

paramfile = 'nemuroParamSets.mat';
A = load(paramfile);

% Set up defaults

params = fieldnames(A.NemParam);
params = params(1:end-1);  % All but setname, at end
setnames = {A.NemParam.setname};
np = length(params);

if nargin == 1 && ischar(varargin{1}) && ismember(varargin{1}, setnames)
    defaultset = varargin{1};
    pv = cell(0);
else
    if nargin > 0 && ischar(varargin{1}) && ismember(varargin{1}, setnames)
        defaultset = varargin{1};
        pv = varargin(2:end);
    else
        defaultset = 'NEMURO Version 1.f90';
        pv = varargin;
    end
end
    
defidx = find(strcmp(setnames, defaultset));

% Parse inputs

p = inputParser;

for ip = 1:np
    p.addParamValue(params{ip}, A.NemParam(defidx).(params{ip}));
end


p.StructExpand = true;
p.KeepUnmatched = true;

p.parse(pv{:});

A = p.Results;

B = p.Unmatched;
