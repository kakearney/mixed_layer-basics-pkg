function varargout = tracerforced(action, varargin)
%TRACERFORCED Generic tracer with deep-water concentration held constant
%
% This module adds a single tracer variable to the mixed_layer model. 
% The tracer has no sources or sinks, but is forced during mixing so that
% the bottom grid cell maintains the same tracer concentration as was
% supplied in the initial profile.
%
% See biomodule.m for function syntax descriptions.  The following 
% fields must be present in the In structure (passed to mixed_layer as
% parameter/value pairs):
%
%   tracer: n x 2 depth profile of tracer, where column 1 gives the depth
%           values (negative down) and column 2 holds the initial
%           concentrations of a generic tracer (mmol m^-3) 

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

bio = interp1(In.tracer(:,1), In.tracer(:,2), Grd.z);
ismixed = true;
Biovars = struct;
bottomval = bio(end);
names = {'tracer', 'a tracer variable', 'no units'};
diagnames = cell(0);


%**************************************************************************

function [newbio, diag] = sourcesink(oldbio, meanqi, temperature, z, dz, ...
                             Biovars, t, dt);
                                                  
newbio = oldbio;
diag = [];

%**************************************************************************

function wsink = vertmove(oldbio, meanqi, temperature, z, dz, Biovars, t, dt)

wsink = zeros(size(oldbio));
