function Con = bioconservecheck(Con, bio, it, id)
%
% Con = bioconservecheck(outfile) 
% Con = bioconservecheck(Con, bio, dz, it, id)
%
% id: 1 = pre-mix
%     2 = post-mix,pre-vert
%     3 = post-vert, pre-source/sink
%     4 = post-all

if nargin == 1
    outfile = Con; % reassignment just to avoid warning
    Con = struct;
    Con.file = regexprep(outfile, '.nc', '_conservecheck.dat');
else

    tot = sum(bsxfun(@times, bio, Con.dz));

    ntot = sum(tot(1:27));          % mol m^-2
    stot = sum(tot([28:29 31]));    % mol m^-2
    ftot = sum(tot([30 32:33]));    % umol m^-2

    if it == 1
        Con.ntot0 = ntot;
        Con.fid = fopen(Con.file, 'wt');
    end
   
    solver = getappdata(0, 'solvercheck');
    if solver == 4
        error('WCE:notConserved', 'Implicit solver');
    end
    
    fprintf(Con.fid, '%f %f %f %f %f %f\n', it, id, ntot, stot, ftot, solver);
    
    if abs(ntot - Con.ntot0) > 0.01
        fclose(Con.fid);
        error('WCE:notConserved', 'N not conserved');
    end
end

