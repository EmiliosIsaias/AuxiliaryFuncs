function [subFolds] = getSubFolds(folder)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
subFolds = dir(folder);
pointFlag = arrayfun(@(x) any(strcmpi(x.name, {'.','..'})), subFolds);
foldFlag = reshape([subFolds.isdir], [], 1);
subFolds(pointFlag | ~foldFlag) = [];
end