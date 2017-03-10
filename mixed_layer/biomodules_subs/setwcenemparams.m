function [Biovars, Np] = setwcenemparams(BioIn, nemidx, nbsv, Grd, Biovars)
%SETWCENEMPARAMS Setup of parameters for mixed_layer wce/nemurokak modules

%-----------------------
% Shared parameters
%-----------------------

% ODE solver

if ischar(BioIn.odesolver)
    Biovars.odesolver = {BioIn.odesolver};
else
    Biovars.odesolver = BioIn.odesolver;
end

% Extended NEMURO parameters with 1:1 correspondence to nemurokak
% (only defined for original-11 NEMURO-derived variables, so reassign to
% proper indices) 

Np = nemparams2arrays(BioIn.NemParam);

Biovars.alpha1                  = Np.alpha1;              % m^-1  
Biovars.alpha2                  = Np.alpha2/1000;         % (m^3/molN)/m
Biovars.usesteele               = Np.usesteele;           % logical
Biovars.RSiN                    = Np.RSiN;                % molSi/molN

Biovars.Iopt                    = zeros(nbsv,1);
Biovars.Iopt(nemidx(1:11))      = Np.Iopt / 0.001433;     % W/m^2   
Biovars.alpha                   = zeros(nbsv,1);
Biovars.alpha(nemidx(1:11))     = Np.alpha * 0.001433;    % (W/m^2)^-1
Biovars.Vmax                    = zeros(nbsv,1);
Biovars.Vmax(nemidx(1:11))      = Np.Vmax;                % s^-1
Biovars.Kno3                    = zeros(nbsv,1);
Biovars.Kno3(nemidx(1:11))      = Np.Kno3 * 1000;         % molN/m^3
Biovars.pusai                   = zeros(nbsv,1);
Biovars.pusai(nemidx(1:11))     = Np.pusai / 1000;        % m^3/molN
Biovars.Knh4                    = zeros(nbsv,1);
Biovars.Knh4(nemidx(1:11))      = Np.Knh4 * 1000;         % molN/m^3 
Biovars.Ksi                     = zeros(nbsv,1);
Biovars.Ksi(nemidx(1:11))       = Np.Ksi * 1000;          % molSi/m^3
Biovars.Kgpp                    = zeros(nbsv,1);
Biovars.Kgpp(nemidx(1:11))      = Np.Kgpp;                % degC^-1
Biovars.gamma                   = zeros(nbsv,1);
Biovars.gamma(nemidx(1:11))     = Np.gamma;               % no unit
Biovars.res0                    = zeros(nbsv,1);
Biovars.res0(nemidx(1:11))      = Np.res0;                % s^-1
Biovars.Kres                    = zeros(nbsv,1);
Biovars.Kres(nemidx(1:11))      = Np.Kres;                % degC^-1
Biovars.Kdec                    = zeros(nbsv,nbsv);
Biovars.Kdec(nemidx(1:11),nemidx(1:11)) = Np.Kdec;        % degC^-1
Biovars.vdec                    = zeros(nbsv,nbsv);
Biovars.vdec(nemidx(1:11),nemidx(1:11)) = Np.vdec;        % s^-1

% Iron parameters: Most scalars, a few assigned to phytoplankton groups
% only

Biovars.Kfe                     = zeros(nbsv,1);
Biovars.Kfe(nemidx(1:2))        = BioIn.Kfe;        % molFe/m^3
Biovars.kfe2n                   = zeros(nbsv,1);
Biovars.kfe2n(nemidx(1:2))      = BioIn.kfe2n;      % molFe/molN
Biovars.fe2nmax                 = zeros(nbsv,1);
Biovars.fe2nmax(nemidx(1:2))    = BioIn.fe2nmax;    % molFe/molN
Biovars.fe2nupfac               = BioIn.fe2nupfac;  % molFe/molN
Biovars.ligbkg                  = BioIn.ligbkg;     % mol/m^3
Biovars.alphascav               = BioIn.alphascav;  % s^-1
Biovars.remineff                = BioIn.remineff;   % no unit
Biovars.kliglo                  = BioIn.kliglo;     % mol lig^-1 kg^-1
Biovars.klighi                  = BioIn.klighi;     % mol lig^-1 kg^-1
Biovars.iofescav                = BioIn.iofescav;   % no unit
Biovars.scavfac                 = BioIn.scavfac;    % m^3/mol Fe
Biovars.scavthresh              = BioIn.scavthresh; % mol Fe m^-3

% PON, Opal, and POFe sink

% if BioIn.isnem
%     Biovars.sink = zeros(1,nbsv);
%     Biovars.sink(1,Biovars.idx.pon)  = -Np.settle(Biovars.idx.pon); 
%     Biovars.sink(1,Biovars.idx.opal) = -Np.settle(Biovars.idx.opal);
%     Biovars.sink(1,Biovars.idx.pofe) = -Np.settle(Biovars.idx.pon); % Same as PON
% else
Biovars.settle = zeros(Grd.nz,nbsv);
Biovars.settle(:,nemidx(1:11)) = repmat(-Np.settle', Grd.nz, 1);        % m/s
if isnan(BioIn.pofesink)
    Biovars.settle(:,Biovars.idx.pofe) = Biovars.settle(:,Biovars.idx.pon); % Same as PON
else
    Biovars.settle(:,Biovars.idx.pofe) = -BioIn.pofesink; 
end


% Prey visibility by predators

zvis = BioIn.preyvis(:,1);
[nr,nc] = size(BioIn.preyvis);

vis = ones(nr,nbsv);
vis(:,1:nc-1) = BioIn.preyvis(:,2:end);
Biovars.preyvis = interp1(zvis, vis, Grd.z);

if any(isnan(Biovars.preyvis(:)))
    error('Visibility spectrum must extend entire water column');
end

% TODO: may need to add in a check to make sure all nekton see the surface
% layer when preying on other nekton, but I'll put that on my
% todo-when-I-need-it list
    

% %-----------------------
% % Unshared parameters
% %-----------------------
% 
% if BioIn.isnem 
%       
%     %-----------------------
%     % Nemurokak-only
%     %-----------------------
%     
% 
% 
% 
% else  
%     
%     %-----------------------
%     % WCE-only
%     %-----------------------
% 
%     % Some variables needed
%     
%     Ewein = ecopathinputcheck(BioIn.Ewein, true);
%     Ep = ecopathlite(Ewein);
%     
%     isabovemld = Grd.z >= BioIn.mld;
%     nlayer = sum(isabovemld);
%     mld = -Grd.zp(nlayer+1);
%     
% 
% end


%-----------------------
% Diapause
%-----------------------

% if Biovars.diapause
%     
% 
%     
% end
