function CheckTargets(MRS_struct)

switch num2str(length(MRS_struct.p.target))
    case '1'
        if any(strcmp(MRS_struct.p.target,{'GABA','Glx','GABAGlx','GSH','Lac','EtOH'}))
            if MRS_struct.p.phantom && strcmp(MRS_struct.p.target,'GABAGlx')
                PassError;
            end
            if MRS_struct.p.HERMES
                error('MRS_struct.p.HERMES is set to 1 in GannetPreInitialise.m. Add a second target metabolite or set flag to 0');
            end
        else
            PassError;
        end
    case '2'
        if any([all(strcmp(MRS_struct.p.target,{'GABAGlx','GSH'})) ...
                all(strcmp(MRS_struct.p.target,{'GABA','GSH'})) ...
                all(strcmp(MRS_struct.p.target,{'Lac','GSH'}))])
            if MRS_struct.p.phantom && any(strcmp(MRS_struct.p.target,'GABAGlx'))
                PassError;
            end
            if ~MRS_struct.p.HERMES
                error('Two target metabolites detected. MRS_struct.p.HERMES must be set to 1 in GannetPreInitialise.m.');
            end
        else
            PassError
        end
    case '3'
        if all(strcmp(MRS_struct.p.target,{'EtOH','GABA','GSH'}))
            if ~MRS_struct.p.HERMES
                error('Three target metabolites detected. MRS_struct.p.HERMES must be set to 1 in GannetPreInitialise.m.');
            end
        else
            PassError;
        end
end

end


function PassError

error('Incorrect targets entered in GannetPreInitialise.m. Check spelling/allowable options.');

end



