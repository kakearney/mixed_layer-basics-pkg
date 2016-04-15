function varargout = tracer(action, varargin)
%TRACER Generic tracer
%
% This module adds a simple tracer variable to the mixed_layer model.  The
% tracer has no sources or sinks.
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
        narginchk(nin(1), nin(1));       
        [out{:}] = init(varargin{:});
        
    case 'sourcesink'
        
        out = cell(1, nout(2));
        narginchk(nin(2), nin(2));       
        [out{:}] = sourcesink(varargin{:});
        
    case 'vertmove'
        
        out = cell(1, nout(3));
        narginchk(nin(3), nin(3));       
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
bottomval = NaN;
names = {'tracer', 'a tracer variable', 'no units'};
diagnames = cell(0);

%**************************************************************************

function [newbio, diag] = sourcesink(oldbio, P, B, G)
                         
newbio = oldbio;
diag = [];

%**************************************************************************

function wsink = vertmove(oldbio, P, B, G)

wsink = zeros(size(oldbio));