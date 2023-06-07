function [outTables] = joinRNTables(inTables)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
for ct = 1:numel(inTables)
    try
        inTables{ct} = cat(1, inTables{ct}{:});
    catch
        for ctt = 1:numel(inTables{ct})
            if ~isstring(inTables{ct}{ctt}.Conditions)
                t2t = cellfun(@(t) any(contains(t, "Delay"),"all"), ...
                    inTables{ct}{ctt}.Conditions);
                inTables{ct}{ctt} = table( ...
                    inTables{ct}{ctt}{t2t,"Conditions"}{:}', ...
                    inTables{ct}{ctt}{t2t,"BehaviourIndices"}{:}', ...
                    'VariableNames', {'Conditions', 'BehaviourIndices'});
            end
        end
        inTables{ct} = cat(1, inTables{ct}{:});
    end
end
outTables = cat(1, inTables{:});
end