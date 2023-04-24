function [foundDir] = recursiveFolderSearch(seedPath, targetFolder)
%RECURSIVEFOLDERSEARCH looks for the given folder (targetFolder) in all
%subfolders from the seed folder (seedPath) and returns the full path of
%the searched folder. It returns the first folder found in the subfolders
%(foundDir). If no folder is found, the function returns an empty variable.

% Emilio Isaias-Camacho @ GrohLab 2023
%%
foundDir = [];
childDirs = dir(seedPath);
expandName = @(x) fullfile(x.folder,x.name);
% Remove the dot folders and files from the search.
pointFlag = arrayfun(@(x) any(strcmpi(x.name, {'.','..'})), childDirs);
fileFlag = reshape([childDirs.isdir], [], 1);
childDirs(pointFlag | ~fileFlag) = [];
childNames = arrayfun(@(x) string(x.name), childDirs);
foundDirFlag = strcmp(childNames, targetFolder);
cf = 0; Ncf = numel(childDirs);
if all(~foundDirFlag)
    % Didn't find the folder. Keep looking in subfolders
    while Ncf > 0 && cf < Ncf && isempty(foundDir)
        foundDir = recursiveFolderSearch(expandName(childDirs(cf+1)), targetFolder);
        cf = cf + 1;
    end
else
    % Found it! Return the path from all instances.
    foundDir = string(expandName(childDirs(foundDirFlag)));
end
end