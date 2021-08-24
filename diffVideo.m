function diffSignal = diffVideo(videoPath)
vidObj = VideoReader(videoPath);
mem = memory;
maxMem = mem.MaxPossibleArrayBytes;
frameByte = vidObj.Height * vidObj.Width * 3;
framesInMem = round((maxMem/frameByte) * 0.85);

end