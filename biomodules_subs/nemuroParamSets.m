function NemParam = nemuroParamSets(setname)
%NEMUROPARAMSETS Reads NEMURO parameter set(s) from file
%
% N = nemuroParamSets
% N = nemuroParamSets(setname)
%
% This function reads in and returns sets of NEMURO input parameters,
% gathered from a variety of sources and converted to a standardized(ish)
% set of units.
%
% Note: this replaces the old nemuroinputparser function.
%
% Input variables:
%
%   setname:    string or cell array of strings, corresponding to available
%               parameter sets.  If included, only the specified parameter
%               sets will be returned; otherwise, all available ones will
%               be returned.  Strings may be:
%
%               'Eslinger et al. Simulation Parameter': from the workshop
%                   report where early versions of NEMURO were tested
%                   (Eslinger DL, Coastal N, Centre S, Hobson S,Werner FE,
%                   Carolina N Final Report of the International Workshop
%                   to Develop a Prototype Lower Trophic Level.
%                   North:1?77).  This writeup represented a work in
%                   progress, and was ambiguous in some places.
%
%               'Eslinger et al. Station P': same citation as previous
%
%               'Eslinger et al. A7': same citation as previous
%
%               'Eslinger et al. Bering' same citation as previous
%
%               'NEMURO Version 1.f90': values taken from the original
%                   NEMURO fortran code
%
%               'Fujii et al. Papa': from Fujii M, Yamanaka Y, Nojiri Y,
%                   Kishi MJ, Chai F (2007) Comparison of seasonal
%                   characteristics in biogeochemistry among the subarctic
%                   North Pacific stations described with a NEMURO-based
%                   marine ecosystem model. Ecol Modell 202:52?67.
%                   Tailored to Ocean Station Papa (50N, 145W)      
%
%               'Fujii et al. A7': same citation as above, tailored to
%                   station A7 (41.5N, 145.5E) 
%
%               'Fujii et al. KNOT': same citation as above, tailored to
%                   station KNOT (44N, 155E)
%
%               'Kishi et al. A7': from Kishi M, Kashiwai M, Ware D, Megrey
%                   B, Eslinger D, Werner F, Noguchiaita M, Azumaya T,
%                   Fujii M, Hashimoto S (2007) NEMURO- a lower trophic
%                   level model for the North Pacific marine ecosystem.
%                   Ecol Modell 202:12?25.
%
%
% Output variables:
%
%   N:          nset x 1 structure array, where each element of the
%               structure array corresponds to a parameter set from a
%               different source. Fields are as follows:  
%
%               alpha1:       Light Dissipation coefficient of sea water  [/m]
%               alpha2:       PS+PL Selfshading coefficientS+PL           [l/molN/m]
%               IoptS:        PS Optimum Light Intensity  S               [ly/min]
%               IoptL:        PL Optimum Light Intensity                  [ly/min]
%               LLN:          Number of sublayer for calculating of Lfc
%               VmaxS:        PS Maximum Photosynthetic rate @0degC       [/s]
%               KNO3S:        PS Half satuation constant for Nitrate      [molN/l]
%               KNH4S:        PS Half satuation constant for Ammonium     [molN/l]
%               PusaiS:       PS Ammonium Inhibition Coefficient          [l/molN]
%               KGppS:        PS Temp. Coeff. for Photosynthetic Rate     [/degC]
%               MorPS0:       PS Mortality Rate @0degC                    [/s]
%               KMorPS:       PS Temp. Coeff. for Mortality               [/degC]
%               ResPS0:       PS Respiration Rate at @0degC               [/s]
%               KResPS:       PS Temp. Coeff. for Respiration             [/degC]
%               GammaS:       PS Ratio of Extracell. Excret. to Photo.    [(nodim)]
%               VmaxL:        PL Maximum Photosynthetic rate @0degC       [/s]
%               KNO3L:        PL Half satuation constant for Nitrate      [molN/l]
%               KNH4L:        PL Half satuation constant for Ammonium     [molN/l]
%               KSiL:         PL Half satuation constant for Silicate     [molSi/l]
%               PusaiL:       PL Ammonium Inhibition Coefficient          [l/molN]
%               KGppL:        PL Temp. Coeff. for Photosynthetic Rate     [/degC]
%               MorPL0:       PL Mortality Rate @0degC                    [/s]
%               KMorPL:       PL Temp. Coeff. for Mortality               [/degC]
%               ResPL0:       PL Respiration Rate at @0degC               [/s]
%               KResPL:       PL Temp. Coeff. for Respiration             [/degC]
%               GammaL:       PL Ratio of Extracell. Excret. to Photo.    [(nodim)]
%               GRmaxS:       ZS Maximum Rate of Grazing PS @0degC        [/s]
%               KGraS:        ZS Temp. Coeff. for Grazing                 [/degC]
%               LamS:         ZS Ivlev constant                           [l/molN]
%               PS2ZSstar:    ZS Threshold Value for Grazing PS           [molN/l]
%               AlphaZS:      ZS Assimilation Efficiency                  [(nodim)]
%               BetaZS:       ZS Growth Efficiency                        [(nodim)]
%               MorZS0:       ZS Mortality Rate @0degC                    [/s]
%               KMorZS:       ZS Temp. Coeff. for Mortality               [/degC]
%               GRmaxLps:     ZL Maximum Rate of Grazing PS @0degC        [/s]
%               GRmaxLpl:     ZL Maximum Rate of Grazing PL @0degC        [/s]
%               GRmaxLzs:     ZL Maximum Rate of Grazing ZS @0degC        [/s]
%               KGraL:        ZL Temp. Coeff. for Grazing                 [/degC]
%               LamL:         ZL Ivlev constant                           [l/molN]
%               PS2ZLstar:    ZL Threshold Value for Grazing PS           [molN/l]
%               PL2ZLstar:    ZL Threshold Value for Grazing PL           [molN/l]
%               ZS2ZLstar:    ZL Threshold Value for Grazing ZS           [molN/l]
%               AlphaZL:      ZL Assimilation Efficiency                  [(nodim)]
%               BetaZL:       ZL Growth Efficiency                        [(nodim)]
%               MorZL0:       ZL Mortality Rate @0degC                    [/s]
%               KMorZL:       ZL Temp. Coeff. for Mortality               [/degC]
%               GRmaxPpl:     ZP Maximum rate of grazing PL @0degC        [/s]
%               GRmaxPzs:     ZP Maximum rate of grazing ZS @0degC        [/s]
%               GRmaxPzl:     ZP Maximum rate of grazing ZL @0degC        [/s]
%               KGraP:        ZP Temp. Coeff. for grazing                 [/degC]
%               LamP:         ZP Ivlev constant                           [l/molN]
%               PL2ZPstar:    ZP Threshold Value for Grazing PL           [molN/l]
%               ZS2ZPstar:    ZP Threshold Value for Grazing ZS           [molN/l]
%               ZL2ZPstar:    ZP Threshold Value for Grazing ZL           [molN/l]
%               PusaiPL:      ZP Preference Coeff. for PL                 [l/molN]
%               PusaiZS:      ZP Preference Coeff. for ZS                 [l/molN]
%               AlphaZP:      ZP Assimilation Efficiency                  [(nodim)]
%               BetaZP:       ZP Growth Efficiency                        [(nodim)]
%               MorZP0:       ZP Mortality Rate @0degC                    [/s]
%               KMorZP:       ZP Temp. Coeff. for Mortality               [/degC]
%               Nit0:         NH4 Nitrification Rate @0degC               [/s]
%               KNit:         NH4 Temp. coefficient for Nitrification     [/degC]
%               VP2N0:        PON eecomp. Rate to Ammonium @0degC         [/s]
%               KP2N:         PON Temp. Coeff. for Decomp. to Ammon.      [/degC]
%               VP2D0:        PON Decomp. Rate to DON @0degC              [/s]
%               KP2D:         PON Temp. Coeff. for Decomp. to DON         [/degC]
%               VD2N0:        DON Decomp. Rate to Ammonium @0degC         [/s]
%               KD2N:         DON Temp. Coeff. for Decomp. to Ammon.      [/degC]
%               VO2S0:        Opal Decomp. Rate to Silicate @0degC        [/s]
%               KO2S:         Opal Temp. Coeff. for Decomp.to Silicate    [/degC]
%               RSiN:         Si/N ratio                                  [molSi/molN]
%               RCN:          C/N ratio                                   [molC/molN]
%               setVPON:      Settling velocity of PON                    [m/s]
%               setVOpal:     Settling velocity of Opal                   [m/s]
%               TNO3d:        Nitrate Concentraion in the Deep Layer      [molN/l]
%               TSiOH4d:      Silicate Concentraion in the Deep Layer     [molSi/l]
%               usesteele:    logical scalar indicating whether to use Steele (1) or
%                             Platt (0) light curves

% Copyright 2017 Kelly Kearney


% Read and parse the data from the Excel spreadsheet

S = warning('off', 'MATLAB:table:ModifiedVarnames');
Ptable = readtable('nemuroParamSets.xls', 'readrownames', true, 'readvariablenames', true);
warning(S);

% Extract the parameter sets into arrays

sets = Ptable{1,3:2:end}; % names of each set
unit = Ptable{3:end,1};

nemparam = Ptable.Properties.RowNames(3:end);

setdata = table2array(Ptable(3:end,2:2:end));
setunit = table2array(Ptable(3:end,3:2:end));

% Filter sets if input provided

if nargin > 0
    [tf,loc] = ismember(setname, sets);
    if ~all(tf)
        str = sprintf('%s,', setname{~tf});
        error('Specified set name(s) (%s) do not match available ones', str(1:end-1));
    end
    setdata = setdata(:,loc);
    setunit = setunit(:,loc);
    sets = sets(loc);
end

[np, ns] = size(setdata);

% First convert the @10degC rates to @0deg rates

for is = 1:ns
    for ip = 1:np
        if regexp(setunit{ip,is}, 'at 10')
            q10 = exp(10*setdata(ip+1,is));
            setdata(ip,is) = setdata(ip,is)/q10;
            setunit{ip,is} = regexprep(setunit{ip,is}, ' at 10', '');
        end
    end
end

% Now convert everything else

unitconvert = cell(0,2);
for is = 1:ns
    unitmatch = [setunit(:,is) unit];
    
    temp = [unitconvert; unitmatch];
    unqtemp = unique(temp);
    [tf,loc] = ismember(temp, unqtemp);
    unq = unique(loc, 'rows');
    unitconvert = unqtemp(unq);
    
end
isbad = any(cellfun(@isempty, unitconvert), 2);
unitconvert = unitconvert(~isbad,:);
issame = strcmp(unitconvert(:,1), unitconvert(:,2));
unitconvert = unitconvert(~issame,:);

convertfac = {...
    '/cm'              '/m'             100    
    '/day'             '/s'             1/86400         
    'W/m^2'            'ly/min'         0.001433     
    'l/molN/cm'        'l/molN/m'       100   
    'l/molN/day'       'l/molN/s'       1/86400
    'l/umol/m'         'l/molN/m'       1e6
    'l/umolN'          'l/molN'         1e6     
    'l/umolN/day'      'l/molN/s'       1e6/86400
    'm/day'            'm/s'            1/86400
    'm^2/W/day'        '/(ly/min)/s'    1/(86400*0.001433)
    'm^2/mmolN'        'l/molN/m'       1e6
    'm^3/mmolN'        'l/molN'         1e6
    'm^3/mmolN/day'    'l/molN/s'       1e6/86400
    'mmolN/m^3'        'molN/l'         1e-6
    'mmolSi/m^3'       'molSi/l'        1e-6
    'umolN/l'          'molN/l'         1e-6
    'umolSi/l'         'molSi/l'        1e-6};

catcell = @(a,b) cellfun(@(x,y) [x ' to ' y], a, b, 'uni', 0);
for is = 1:ns
    temp1 = catcell(setunit(:,is), unit);
    temp2 = catcell(convertfac(:,1), convertfac(:,2));
    [tf,icf] = ismember(temp1, temp2);
    setdata(tf,is) = setdata(tf,is) .* cell2mat(convertfac(icf(tf),3));
end

% Check whether uses Steele or Platt light curve

[tf, loc] = ismember({'IoptS'; 'IoptL'; 'alphaS'; 'alphaL'}, nemparam);
lightParams = setdata(loc,:)';

issteele = ~any(isnan(lightParams(:,1:2)),2) & all(isnan(lightParams(:,3:4)),2);
isplatt = all(isnan(lightParams(:,1:2)),2) & ~any(isnan(lightParams(:,3:4)),2);
if ~all(issteele | isplatt)
    error('At least one model doesn''t have correct light parameters');
end

nemparam = [nemparam; {'usesteele'}];
setdata = [setdata; issteele'];

nemparam = [nemparam; {'setname'}];
setdata = [num2cell(setdata); sets];

NemParam = cell2struct(setdata, nemparam, 1);

       