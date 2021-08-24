function diffSignal = diffVideo(videoPath)
vidObj = VideoReader(videoPath);
if vidObj.hasFrame
    frame = getFrame(vidObj);
else
    fprintf(2, "The given path contains no frames. Maybe it wasn't a video/n")
    return
end
diffSignal = zeros(vidObj.NumFrames - 1,1); vidCount = 1;
while vidObj.hasFrame
    nextFrame = getFrame(vidObj);
    frameDiff = diff(cat(3, frame, nextFrame), 1, 3);
    % Log10 of the squared sum from all pixels
    diffSignal(vidCount) = log10(sum(frameDiff(:).^2));
    frame = nextFrame; vidCount = vidCount + 1;
end
end

function frame = getFrame(vidObj)
frame = vidObj.readFrame; frame = mean(frame, 3, 'omitnan');
end