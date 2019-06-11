function CheckTargets(MRS_struct)

switch num2str(length(MRS_struct.p.target))
    case '1'
        if any(strcmp(MRS_struct.p.target,{'GABA','Glx','GABAGlx','GSH','Lac','EtOH'}))
            % do nothing
        else
            PassError;
        end
    case '2'
        if all(strcmp(MRS_struct.p.target,{'GABAGlx','GSH'})) || all(strcmp(MRS_struct.p.target,{'Lac','GSH'}))
            % do nothing
        else
            PassError
        end
    case '3'
        if all(strcmp(MRS_struct.p.target,{'EtOH','GABA','GSH'}))
            % do nothing
        else
            PassError;
        end
end

end


function PassError

error('Incorrect targets entered in GannetPreInitialise.m. Check order and/or spelling.');

end



