function G = ecopath2graph(Ewein)
%ECOPATH2GRAPH Create graph objects for Ecopath model
%
% [G,H] = ecopath2graph(Ewein)
%
% Input variables:
%
%   Ewein: ecopath input structure
%
% Output variables:
%
%   G:      graph object describing all fluxes within Ecopath model.

% Copyright 2015 Kelly Kearney

%--------------------------
% Create basic graph object
%--------------------------

% Main nodes: groups, fleets, and out-of-system sources and sinks 

Ewein = ecopathinputcheck(Ewein, true);
Ep = ecopathlite(Ewein, 'silent', true);

if ~isfield(Ewein, 'name')
    Ewein.name = strtrim(cellstr(num2str((1:Ewein.ngroup)', 'group%d')));
end
if ~isfield(Ewein, 'fleet')
    Ewein.name = strtrim(cellstr(num2str((1:Ewein.ngear)', 'fleet%d')));
end

outofsys = {'Import/export', 'Respiration', 'Landings', 'Primary production'}';
noos = length(outofsys);

nodenames = [Ewein.name; Ewein.fleet; outofsys];
nnode = length(nodenames);

G = digraph(Ep.flow, nodenames);

% Add a few properties from Ecopath calculations
% Type: 0 = consumer, 1 = producer, 2 = detritus, 3 = fleet, 
%       4 = out-of-system

G.Nodes.B = [Ep.b; sum(Ewein.landing,1)'; zeros(noos,1)]; % TODO: assuming /yr
G.Nodes.type = [ceil(Ewein.pp); ones(Ewein.ngear,1)*3; ones(noos,1)*4]; 

% Redo trophic level calculations, including groups and gear

adj = zeros(Ewein.ngroup+Ewein.ngear);
adj(1:Ewein.ngroup,1:Ewein.nlive) = Ep.q0(1:Ewein.ngroup,1:Ewein.nlive);
adj(1:Ewein.ngroup,(1:Ewein.ngear)+Ewein.ngroup) = Ewein.landing+Ewein.discard;
dc = bsxfun(@rdivide, adj, sum(adj,1));    
dc(isnan(dc)) = 0;

pp = [Ewein.pp; zeros(Ewein.ngear,1)];

[pp, isrt] = sort(pp);
dc = dc(isrt,isrt);

tl = trophiclevel(dc, pp, Ewein.nlive+Ewein.ngear, Ewein.ngroup+Ewein.ngear);
[~,loc] = ismember(1:(Ewein.ngroup+Ewein.ngear), isrt);
G.Nodes.TL = [tl(loc); nan(noos,1)];


