function smudMeanPx = diffVideo(videoPath)
if ~exist(videoPath, "file")
    fprintf(1, 'Given path does not exist!\n')
    return
end
dataDir = fileparts(videoPath);
dt = getDates(videoPath,'roller');
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
        possFramesInMem = round((maxMem/frameByte) * 0.73);
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
    save(foName, 'fr', "sPx")
else
    fprintf(1, 'Video movement file exists!\nNo file nor output created\n')
end
