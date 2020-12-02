function saveFigure(figSave,absBaseFilePath,vectorFlag)
if ~exist('vectorFlag','var')
    vectorFlag = 1;
end
absBaseFilePath = string(absBaseFilePath); fileExt = [".fig";".pdf";".emf"];
eidff = arrayfun(@(x) ~exist(absBaseFilePath + x,'file'),...
    fileExt);
if eidff(1)
    savefig(figSave, absBaseFilePath + fileExt(1))
end
if vectorFlag
    figSave = configureFigureToPDF(figSave);
end
if eidff(2)
    print(figSave, absBaseFilePath + fileExt(2),'-dpdf','-fillpage')
end
if eidff(3)
    print(figSave, absBaseFilePath + fileExt(3),'-dmeta')
end