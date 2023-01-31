% Extend file name from file structure
exfn = @(x) fullfile(x.folder, x.name);
getCondNames = @(y) arrayfun(@(x) string(x.ConditionName), y);
findStr = @(x, y) contains(x, y, "IgnoreCase", true);
parentDir = "Z:\Emilio\SuperiorColliculusExperiments\Roller\Batch11_ephys.MC";
summFiles = dir(fullfile(parentDir, "*","*","*","Behaviour/Simple summary.mat"));
for cf = summFiles'
    % summary file name
    sfn = exfn(cf); fprintf(1, "Processing %s...\n", sfn)
    load(sfn, 'summStruct'); condNames = getCondNames(summStruct);
    ctrlSub = find(findStr(condNames, "Control"));
    lsrSub = find(findStr(condNames, "Delay"));
    if isempty(lsrSub)
        fprintf(1,'No laser condition found in this experiment... ')
        fprintf(1, 'Skipping\n')
        continue
    end
    
end