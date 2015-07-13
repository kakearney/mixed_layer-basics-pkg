function Bio = initbiorelax(Grd, Bio, In)
%INITBIORELAX Sets up relaxation for biological variables
%
% Fields added to Bio structure:
% 
%   hasrelax:   nbsv x 1 logical array, true if relaxation data provided
%
%   Relax:      1 x nbsv 1 x 1 structure of data for bio relaxation
%               interpolation.  Fields empty for non-relaxed variables.
%               t:      nt x 1 array, times corresponding to
%                       columns of data (s from sim start time) 
%               z:      nz x 1 array, depths corresponding to the
%                       rows of data (m, neg down) 
%               data:   nz x nt array, biological relaxation profiles
%                       (in units defined bio biological module) 
%
%   

%-------------------------
% Relaxation
%-------------------------

% Determine the names of relaxation input variables (xxxrelax, where xxx
% are the short names of all biological state variables)

relaxvar = cellfun(@(x) [x 'relax'], Bio.names(:,1), 'uni', 0);
inflds = fieldnames(In);

% Indicator variable for relaxation

Bio.hasrelax = ismember(relaxvar, inflds);

% Set up interpolation grids

nb = size(Bio.names, 1);

for ib = 1:nb
    if Bio.hasrelax(ib)
        [Bio.Relax(ib), str] = initinterpdata('both', In.(relaxvar{ib}), Grd);
    
        if ~isempty(str) && In.verbose
            fprintf('  Missing %s relaxation data: %s\n', Bio.names{ib,1}, str);
        end
    end
end


%-------------------------
% Outside flux
%-------------------------

% Determine the names of relaxation input variables (xxxrelax, where xxx
% are the short names of all biological state variables)

extrafluxvar = cellfun(@(x) [x 'flux'], Bio.names(:,1), 'uni', 0);
inflds = fieldnames(In);

% Indicator variable for relaxation

Bio.hasflux = ismember(extrafluxvar, inflds);

% Set up interpolation grids

for ib = 1:nb
    if Bio.hasflux(ib)
        if isempty(In.(extrafluxvar{ib}))
            Bio.hasflux(ib) = false; % B/c I don't have an earlier error check for empty
        else
            [Bio.ExtraFlux(ib), str] = initinterpdata('both', In.(extrafluxvar{ib}), Grd);

            if ~isempty(str) && In.verbose
                fprintf('  Missing %s outside flux data: %s\n', Bio.names{ib,1}, str);
            end
        end
    end
end
