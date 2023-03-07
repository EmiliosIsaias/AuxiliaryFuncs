function rewriteClusterWaveforms(batchDir)
expandName = @(x) fullfile(x.folder, x.name);
    function renameWaveformFile()
        
        wFiles = dir(fullfile(expandName(cef), '*_waveforms.mat'));
        if ~isempty(wFiles)
            movefile(expandName(wFiles), fullfile(wFiles.folder, ...
                insertBefore(wFiles.name,".","BACKUP")));
        end
    end
ephDirs = dir(fullfile(batchDir,"\**\ephys*"));
ephDirs(~[ephDirs.isdir]) = [];
for cef = ephDirs(:)'
    renameWaveformFile()
    ciFile = fullfile(expandName(cef), "cluster_info.tsv");
    try
        clInfo = getClusterInfo(ciFile);
        varNames = clInfo.Properties.VariableNames;
        varFlag = ismember(varNames, {'cluster_id', 'id'});
        gclID = cellstr(clInfo{clInfo.ActiveUnit==1, varNames(varFlag)});
        fprintf(1, "Processing %s... ", expandName(cef))
        getClusterWaveform(gclID, expandName(cef));
        fprintf(1, "Done!\n")
    catch ME
        disp(ME)
        continue
    end
end

end

