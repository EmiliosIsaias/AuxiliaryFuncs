function [newStr] = capitalizeFirst(inputStr)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
newStr = insertBefore(eraseBetween(inputStr,1,1),...
    1,upper(extractBefore(inputStr,2)));
end

