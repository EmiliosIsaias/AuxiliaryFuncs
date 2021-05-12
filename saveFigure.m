function saveFigure(figSave,absBaseFilePath, vectorFlag, ovrWriteFlag)
if ~exist('vectorFlag','var')
    vectorFlag = 1;
end
if ~exist('ovrWriteFlag','var')
    ovrWriteFlag = false;
end

absBaseFilePath = string(absBaseFilePath); fileExt = [".fig";".pdf";".emf"];
eidff = arrayfun(@(x) ~exist(absBaseFilePath + x,'file'),...
    fileExt);
if eidff(1) || ovrWriteFlag
    savefig(figSave, absBaseFilePath + fileExt(1))
end
if vectorFlag
    figSave = configureFigureToPDF(figSave);
end
if eidff(2) || ovrWriteFlag
    print(figSave, absBaseFilePath + fileExt(2),'-dpdf','-fillpage')
end
if eidff(3) || ovrWriteFlag
    print(figSave, absBaseFilePath + fileExt(3),'-dmeta')
end