function saveFigure(figSave,absBaseFilePath, vectorFlag, ovrWriteFlag)
if ~exist('vectorFlag','var')
    vectorFlag = true;
end
if ~exist('ovrWriteFlag','var')
    ovrWriteFlag = false;
    askOvr = true;
else
    askOvr = false;
end

absBaseFilePath = string(absBaseFilePath); fileExt = [".fig";".pdf";".emf"];
eidff = arrayfun(@(x) ~exist(absBaseFilePath + x,'file'),...
    fileExt);
questionOpts = {'Overwrite?','Yes','No','No'};
funCell = {@savefig, @print, @print};
savingOpts = {{},{'-dpdf','-fillpage'},{'-dmeta'}};
for cft = 1:3
    if eidff(cft) || ovrWriteFlag
        if askOvr && ~eidff(cft)
            ovrAns = questdlg(sprintf('Overwrite %s?', absBaseFilePath(cft)),...
                questionOpts{:});
            if strcmpi(ovrAns,'No')
                continue
            end
        end
        if isempty(savingOpts{cft})
            funCell{cft}(figSave, absBaseFilePath + fileExt(cft))
        else
            funCell{cft}(figSave, absBaseFilePath + fileExt(cft),...
                savingOpts{cft}{:})
        end
    end
    if vectorFlag && cft == 1 
        figSave = configureFigureToPDF(figSave);
    end
end