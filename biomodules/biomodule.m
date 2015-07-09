function varargout = biomodule(action, varargin)
%BIOMODULE Template biological module for mixed layer model
%
% [bioinit, ismixed, bottomval, Biovars, names] = ...
%    biomodule('init', In, Grd);
% [newbio, diag] = ...
%    biomodule('sourcesink', oldbio, PhysParams, BioParams, GrdParams);
% wsink = ...
%    biomodule('vertmove', oldbio, PhysParams, BioParams, GrdParams);
%
% This is a template for a biological module for the mixed layer model.
% The module is called in three different modes during the mixed_layer
% simulation: 'init', 'sourcesink', and 'vertmove'.  
% 
% The initialize ('init') mode sets up several different variables
% associated with the biology of the model, including: 
%   - setting initial depth profiles of all biological state variables
%   - determining whether each state variable will be mixed via the
%     diffusive mixing scheme used in the model
%   - setting any deep water forcing of each state variable
%   - adding any additional variables that will be needed to run the 
%     'sourcesink' and/or 'vertmove' modes of the module. 
%   - setting the names of any dignostic variables that will be calculated
%     in the 'sourcesink' mode. 
%
% The source and sink mode ('sourcesink') calculates the change in biomass
% of each state variable over a single time step of the model.  This mode
% receives as input the biomass of all biological state variables, the mean
% solar radiation, and the temperature profile at the current time step, as
% well as the time and depth grid variables.  If you require additional
% variables to calculate changes in biomass, these variables should be
% included in the Biovars structure created during initialization.
%
% The vertical movement ('vertmove') mode calculates speed of vertical
% movement for each state variable.  This can be used to introduce sinking
% or any other vertical movement unrelated to turbulent mixing.
%
% Input variables for 'init' mode:
%
%   In:         1 x 1 structure holding user-supplied info.  
%               
%               See each individual function in this folder for the
%               particular variables required by each module.
%
%   Grd:        1 x 1 structure, see mixed_layer code for fields
%
% Output variables for 'init' mode:
%
%   bioinit:    ndepth x nbsv array, initial depth profiles for all
%               biological state variables
%
%   ismixed:    1 x nbsv logical array, indicates whether each biological
%               state variable undergoes physical mixing.
%
%   bottomval:  1 x nbsv array, value used to force bottom layer during
%               mixing.  NaN indicates no forcing.
%
%   Biovars:    1 x 1 structure holding additional parameters
%
%   names:      1 x nbsv cell array of strings, names for each state
%               variable,  Names must begin with a letter and contain only
%               letters, numbers, and underscores.
%
%   diagnames:  1 x ndiag cell array of strings, names for each diagnostic
%               variable.  Names must begin with a letter and contain only
%               letters, numbers, and underscores.
%
% Input variables for 'sourcesink' and 'vertmove' mode:
%
%
%   oldbio:     ndepth x nbsv array, profiles of biological state variables
%               at current time step
%
%   PhysParams: 1 x 1 structure holding water-column physical properties
%               relevant to biological calculations:
%
%               par:    photosynthetically active radiation (W m^-2)
%
%               par24:  photsyntheticaly active radiation, averaged over
%                       the previous 24 hours (W m^-2)  
%
%               kpar:   attenuation coefficient for PAR (m^-1)
%
%               T:      nz x 1 array of temperature
%
%               S:      nz x 1 array of salinity
%
%               Sig:    nz x 1 array of density
%
%   BioParams:  1 x 1 structure holding all biology-related parameters set
%               up in the 'init' phase.  The fields of this structure will
%               vary based on the specific biological module being used.
%
%   GrdParams:  1 x 1 structure holding parameters related to the
%               spatiotemporal model grid:
%
%               z:      ndepth x 1 array, depth grid for model (m)
%
%               dz:     depth grid cell size (m)
%
%               t:      current time (s)
%
%               dt:     time step (s)
%
% Output variables for 'sourcesink' mode:
%
%   newbio:     ndepth x nbsv array, profiles of biological state variables
%               after stepping forward by one time step
%
%   diag:       ndepth x ndiag array, values of each diagnostic variable at
%               each depth
%
% Output variables for 'vertmove' mode:
%
%   wsink:      ndepth x nbsv array, velocity of vertical movement, where
%               positive values indicate movement towards the surface (m/s)

% Copyright 2008 Kelly Kearney

nin(1) = nargin(@init) + 1;
nin(2) = nargin(@sourcesink) + 1;
nin(3) = nargin(@vertmove) + 1;

nout(1) = nargout(@init);
nout(2) = nargout(@sourcesink);
nout(3) = nargout(@vertmove);

switch action
    case 'init'
        
        out = cell(1, nout(1));
        error(nargchk(nargin,nin(1),nin(1)));        
        [out{:}] = init(varargin{:});
        
    case 'sourcesink'
        
        out = cell(1, nout(2));
        error(nargchk(nargin,nin(2),nin(2)));        
        [out{:}] = sourcesink(varargin{:});
        
    case 'vertmove'
        
        out = cell(1, nout(3));
        error(nargchk(nargin,nin(3),nin(3)));        
        [out{:}] = vertmove(varargin{:});
        
    otherwise
        
        error('Invalid action for biological module');
end

[varargout{1:nargout}] = out{:};
        
%**************************************************************************

function [bio, ismixed, bottomval, Biovars, names, diagnames] = init(In, Grd)

% MODIFY CODE HERE

bio = ones(size(Grd.z));
ismixed = true;
Biovars = struct;
bottomval = NaN;
names = {'testState', 'test state veriable', 'no units'};
diagnames = {'testDiagnostic', 'test diagnostic variable', 'no units'};

%**************************************************************************

function [newbio, diag] = sourcesink(oldbio, P, B, G)

% MODIFY CODE HERE
                         
newbio = oldbio;
diag = rand(size(oldbio,1), 1);

%**************************************************************************

function wsink = vertmove(oldbio, P, B, G)

% MODIFY CODE HERE

wsink = zeros(size(oldbio));