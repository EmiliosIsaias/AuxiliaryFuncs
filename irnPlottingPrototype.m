mcirn = nan(sum(mcirnTables{:,"Conditions"}=="Control Puff"), 3);
rcnt = 0;
for cr = 1:size(mcirnTables,1)
    if mcirnTables{cr,"Conditions"} == "Control Puff"
        rcnt = rcnt + 1;
        ccnt = 1;
    else
        if contains(mcirnTables{cr,"Conditions"}, '+ L')
            ccnt = 3;
        else
            ccnt = 2;
        end
    end
    mcirn(rcnt, ccnt) = mcirnTables{cr,"BehaviourIndices"};
end

%%
for ct = 1:numel(mcgrnTables)
    try
        mcgrnTables{ct} = cat(1, mcgrnTables{ct}{:});
    catch
        for ctt = 1:numel(mcgrnTables{ct})
            if ~isstring(mcgrnTables{ct}{ctt}.Conditions)
                t2t = cellfun(@(t) any(contains(t, "Delay"),"all"), ...
                    mcgrnTables{ct}{ctt}.Conditions);
                mcgrnTables{ct}{ctt} = table( ...
                    mcgrnTables{ct}{ctt}{t2t,"Conditions"}{:}', ...
                    mcgrnTables{ct}{ctt}{t2t,"BehaviourIndices"}{:}', ...
                    'VariableNames', {'Conditions', 'BehaviourIndices'});
            end
        end
        mcgrnTables{ct} = cat(1, mcgrnTables{ct}{:});
    end
end
mcgrnTables = cat(1, mcgrnTables{:});