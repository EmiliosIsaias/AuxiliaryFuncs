expandName = @(x) string(fullfile(x.folder, x.name));
dataDir = "Z:\Emilio\SuperiorColliculusExperiments\Roller\Batch13_beh";
% All mat files, skipping the arduino files - no video dependence
resultFiles = dir(fullfile(dataDir, "*\*\*\*.mat"));
resultFiles(arrayfun(@(s) contains(s.name, 'arduino', ...
    'IgnoreCase', true), resultFiles)) = [];
% Looping through all intensity directory rather then each mat file.
intDirs = unique(arrayfun(@(s) string(s.folder), resultFiles));
for cid = intDirs(:)'
    fidFlag = ~isempty(dir(fullfile(cid, 'FrameID*.csv')));
    if fidFlag
        bckDir = fullfile(cid, 'BACKUP');
        if ~exist(bckDir, 'dir')
            if ~mkdir(bckDir)
                fprintf(1, "Couldn't create BACKUP folder in %s\n", cid)
                continue
            end
        end
        arrayfun(@(crf) ...
            movefile(expandName(crf), fullfile(bckDir, crf.name)), ...
            resultFiles(arrayfun(@(x) ismember(x.folder, cid), resultFiles)))
        behRes = analyseBehaviour(cid, 'verbose', false);
        [defProb, ~, behIdxFig] = createBehaviourIndex(behRes);
        fprintf(1, "DP: %.3f\n", defProb); pause
        close all
    else
        fprintf(1, "No Frame ID in %s\n", cid)
    end
end