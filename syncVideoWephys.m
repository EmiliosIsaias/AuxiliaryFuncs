function smudMeanPx = syncVideoWephys(videoPath, varargin)
%% Input validation
p = inputParser;
addRequired(p, 'videoPath', @(x) exist(x, "dir"))
addOptional(p, 'verbose', false, @(x) islogical(x) & numel(x)==1)

p.parse(videoPath, varargin{:})

videoPath = p.Results.videoPath;
verbose = p.Results.verbose;
%% Auxiliary functions
fnOpts = {'UniformOutput', false};
search4This = @(x) dir(fullfile(videoPath, x));
getFullPath = @(x) string(fullfile(x.folder, x.name));
getNameFromVid = @(x) string(fullfile(x.Path, x.Name));
%% Searching for videos and their time axis
vidPaths = search4This('*.avi');
frmIDPaths = search4This('FrameID*.csv');
trigPaths = search4This('TriggerSignals*.bin');
Nv = numel(vidPaths); Nfi = numel(frmIDPaths); Nt = numel(trigPaths);
if Nv ~= Nfi || Nv ~= Nt
    fprintf(1, "No full correspondance in video and frame ID paths!\n")
    fprintf(1, "Cannot continue.\n")
    return
end
vPaths = arrayfun(getFullPath, vidPaths);
vidObj = arrayfun(@VideoReader, vPaths, fnOpts{:});
fr = cellfun(@(x) x.FrameRate, vidObj);
tPaths = arrayfun(getFullPath, trigPaths);
fiPaths = arrayfun(getFullPath, frmIDPaths);
vidTx = arrayfun(@(x) readtable(x, "Delimiter", ","), fiPaths, fnOpts{:});
vidTx = cellfun(@(x) x.Var2/1e9, vidTx, fnOpts{:});
ephTx = arrayfun(@(x) readTriggerFile(x), tPaths, fnOpts{:});

%% Searching for laser
for cv = 1:Nv
    %dt = getDates(vPaths(cv), 'roller');
    %foName = fullfile(dataDir, "VideoMovements" + string(dt));
    if vidObj{cv}.hasFrame
        if verbose; fprintf(1, "Getting number of frames... "); end
        Nf = vidObj{cv}.NumFrames; 
        if verbose; fprintf(1, "Done!\n"); end
        frameByte = vidObj{cv}.Height * vidObj{cv}.Width * 3;
        frCount = 0; auxFrame = [];
        while vidObj{cv}.hasFrame && frCount < Nf
            mem = memory; maxMem = mem.MemAvailableAllArrays;
            possFramesInMem = floor((maxMem/frameByte) * 0.7);
            if frCount + possFramesInMem > Nf
                possFramesInMem = Nf - frCount;
            end
            fprintf(1, 'Reading %d/%d frames...\n', frCount + possFramesInMem, Nf)
            %frames = vidObj{cv}.read(frCount + [1, possFramesInMem]);
            %frames(:,:,[2,3],:) = [];
            %frames = squeeze(frames);
            %frames = cat(3, auxFrame, frames);
            %auxFrame = frames(:,:,end);
            %frCount = frCount + possFramesInMem;
        end
    else
        fprintf(1, "Video %s has no frames!\n", vidObj{cv}.Name)
        return
    end

end
dataDir = fileparts(videoPath);
dt = getDates(videoPath, 'roller');
foName = fullfile(dataDir, "VideoMovements" + string(dt));
if ~exist(foName,"file")
    vidObj = VideoReader(videoPath); Nf = 0; fr = 0;
    smudMeanPx = [];
    if vidObj.hasFrame
        fprintf(1, 'Calculating total number of frames...\n');
        Nf = vidObj.NumFrames; fr = vidObj.FrameRate;
        smudMeanPx = zeros(Nf, 1, 'single');
    else
        fprintf(2, "The given video has no frames!/n")
        return
    end
    frCount = 0; auxFrame = [];
    while vidObj.hasFrame && frCount < Nf
        mem = memory; maxMem = mem.MemAvailableAllArrays;
        frameByte = vidObj.Height * vidObj.Width * 3;
        possFramesInMem = floor((maxMem/frameByte) * 0.7);
        if frCount + possFramesInMem > Nf
            possFramesInMem = Nf - frCount;
        end
        fprintf(1, 'Reading %d/%d frames...\n', frCount + possFramesInMem, Nf)
        frames = vidObj.read(frCount + [1, possFramesInMem]);
        frames(:,:,[2,3],:) = [];
        frames = squeeze(frames);
        frames = cat(3, auxFrame, frames);
        deltaFrame = diff(frames, 1, 3); meanPx = squeeze(mean(deltaFrame, [1, 2]));
        smudMeanPx(1+frCount: frCount+length(meanPx)) = movmedian(meanPx, 65);
        % Auxiliary variables
        auxFrame = frames(:,:,end);
        frCount = frCount + possFramesInMem;
    end
    sPx = smudMeanPx;
    save(foName, 'fr', 'sPx')
else
    fprintf(1, 'Video movement file exists!\nNo file nor output created\n')
end
end
function [trig] = readTriggerFile(tfName)
fID = fopen(tfName,'r'); trig = fread(fID, [2,Inf], 'uint16=>int32');
[~] = fclose(fID); trig = int16(trig - median(trig, 2));
trig(1,:) = -trig(1,:);
end