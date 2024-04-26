function smudMeanPx = diffVideo(videoPath)

smudMeanPx = zeros(1, 1, 'single');

if ~exist(videoPath, "file")
    fprintf(1, 'Given path does not exist!\n')
    return
end
dataDir = fileparts(videoPath);
dt = getDates(string( videoPath ) ,'roller');
foName = fullfile(dataDir, "VideoMovements" + string(dt));

if ~exist(foName,"file")
    vidObj = VideoReader(videoPath);
    ht = vidObj.Height; wd = vidObj.Width;
    frameByte = ht * wd * 3;
    try
        [~, mem] = memory; maxMem = mem.PhysicalMemory.Available;
        bufferFrames = floor( (maxMem/frameByte) * 0.5 );
    catch
        % No memory module installed. Considering 32 GB of RAM memory
        bufferFrames = floor( (32e9/frameByte) * 0.5 );
    end

    fprintf(1, 'Calculating total number of frames... ');
    Nf = vidObj.NumFrames; fr = vidObj.FrameRate;
    smudMeanPx = zeros( Nf-1, 1, 'single' );
    fprintf(1, '%d\n', Nf )
    
    frCount = 0; auxFrame = [];
    while frCount < Nf

        if frCount + bufferFrames > Nf
            bufferFrames = Nf - frCount;
        end
        fprintf(1, 'Reading %d to %d out of %d frames...\n', ...
            frCount + [1, bufferFrames], Nf)
        frames = read( vidObj, frCount + [1, bufferFrames] );
        frames = cat(4, auxFrame, frames);
        deltaFrame = diff( frames(:, :, 1, :), 1, 4 ); 
        meanPx = squeeze( mean( deltaFrame, [1, 2] ) );
        smudMeanPx( frCount + (1:(bufferFrames-1)) ) = ...
            movmedian( meanPx, 65 );

        % Auxiliary variables
        auxFrame = frames(:, :, :, end);
        frCount = frCount + bufferFrames;
    end
    sPx = smudMeanPx;
    save(foName, 'fr', 'sPx')
else
    fprintf(1, 'Video movement file exists!\nNo file nor output created\n')
end
