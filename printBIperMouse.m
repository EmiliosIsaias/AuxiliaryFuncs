function printBIperMouse(miceStr, cDir)

Nm = numel(miceStr);
Ns = arrayfun(@(m) numel(m.Sessions), miceStr);
[~, fName] = fileparts(cDir);
fPath = fullfile(cDir, fName);
fPath = string(fPath)+".tsv";
if exist(fPath, 'file')
    fprintf(1, 'File exists!\nOverwrite?\n')
    
    return
end
fID = fopen(fPath, "w");
if fID < 3
    fprintf(1, 'Error creating the file!\n')
    return
end
fprintf(fID, "Mouse\tInfo\tDate\t#Session\tConditions\n");
try
    for cm = 1:Nm
        for cs = 1:Ns(cm)
            tmpTbl = miceStr(cm).Sessions(cs).DataTable;
            if strcmpi(miceStr(cm).Sessions(cs).Type, 'multi')
                puffTbl = cellfun(@(x) any(contains(x, 'Puff')), ...
                    miceStr(cm).Sessions(cs).DataTable.Conditions);
                tmpTbl = table(miceStr(cm).Sessions(cs).DataTable{puffTbl,1}{:}(:),...
                    miceStr(cm).Sessions(cs).DataTable{puffTbl,2}{:}(:), ...
                    'VariableNames', {'Conditions', 'BehaviourIndex'});
            end
            Nc = size(tmpTbl,1);
            condBIs = "";
            for cc = 1:Nc-1
                condBIs = condBIs + sprintf("%s\t%.3f\t", ...
                    tmpTbl{cc,:});
            end
            condBIs = condBIs + sprintf("%s\t%.3f\n", ...
                tmpTbl{Nc,:});
            % Writing mouse, 'structure', session date, session number, and conditions
            fprintf(fID, "%s\t%s\t%s\t%d\t%s", miceStr(cm).Name, ...
                miceStr(cm).Structure, miceStr(cm).Sessions(cs).Date, cs,...
                condBIs);
        end
    end
    fclose(fID);
catch
    fclose(fID);
end