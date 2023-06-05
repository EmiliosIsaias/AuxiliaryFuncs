function [foundDir] = recursiveFolderSearch(seedPath, targetString, varargin)
%RECURSIVEFOLDERSEARCH looks for the given folder (targetFolder) in all
%subfolders from the seed folder (seedPath) and returns the full path of
%the searched folder. It returns the first folder found in the subfolders
%(foundDir). If no folder is found, the function returns an empty variable.

% Emilio Isaias-Camacho @ GrohLab 2023
%%
p = inputParser;
istxt = @(x) ischar(x) | isstring(x);
isfolder = @(x) istxt(x) & exist(x, 'dir');
addRequired(p, 'seedPath', isfolder)
addRequired(p, 'targetString', istxt)
addParameter(p, 'SearchType', 'string', istxt)

parse(p, seedPath, targetString, varargin{:})

seedPath = p.Results.seedPath;
targetString = p.Results.targetString;
sType = p.Results.SearchType;

foundDir = [];
childDirs = dir(seedPath);
expandName = @(x) fullfile(x.folder,x.name);
% Remove the dot folders and files from the search.
pointFlag = arrayfun(@(x) any(strcmpi(x.name, {'.','..'})), childDirs);
fileFlag = reshape([childDirs.isdir], [], 1);
childDirs(pointFlag | ~fileFlag) = [];
childNames = arrayfun(@(x) string(x.name), childDirs);
if ~isempty(childNames)
    switch sType
        case 'string'
            foundDirFlag = strcmp(childNames, targetString);
        case 'expression'
            foundDirFlag = cellfun(@(x) ~isempty(x), ...
                (regexp(childNames, targetString, 'match')));
    end
    cf = 0; Ncf = numel(childDirs);
    if all(~foundDirFlag)
        % Didn't find the folder. Keep looking in subfolders
        while Ncf > 0 && cf < Ncf && isempty(foundDir)
            foundDir = recursiveFolderSearch(expandName(childDirs(cf+1)), ...
                targetString, 'SearchType', sType);
            cf = cf + 1;
        end
    else
        % Found it! Return the path from all instances.
        foundDir = arrayfun(@(x) string(expandName(x)), childDirs(foundDirFlag));
    end
end
end