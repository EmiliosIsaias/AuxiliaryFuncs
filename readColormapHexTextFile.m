function [clrMap] = readColormapHexTextFile(filepath)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
auxSubs = 1:2:5;
try
    fInfo = dir(filepath);
    % Weird that the file system omits 2 bytes from the file size. Hence
    % the '+2'. Then divided by 8 because that is the number of characters
    % in a line.
    Nclrs = (fInfo.bytes + 2)/8;
    fID = fopen(filepath, "r");
    if fID >= 3
        cnt = 1; clrMap = zeros(ceil(Nclrs), 3);
        while ~feof(fID)
            ln = fgetl(fID);
            clrMap(cnt,:) = arrayfun(@(x,y) hex2dec(ln(x:y)), auxSubs,...
                auxSubs+1); cnt = cnt + 1;
        end
        [~] = fclose(fID);
    else
        fprintf(1, 'Unable to open file: %s\n', filepath)
        fprintf(1, 'No colormap created...\n')
    end
catch
    fprintf(1, '%s not recognizable as a filepath!\n', filepath)
    fprintf(1, 'Unable to get a colormap...\n')
end
end

