function saveFigure(figSave,absBaseFilePath)
absBaseFilePath = string(absBaseFilePath); fileExt = [".fig";".pdf";".emf"];
eidff = arrayfun(@(x) ~exist(absBaseFilePath + x,'file'),...
    fileExt);
if eidff(1)
    savefig(figSave, absBaseFilePath + fileExt(1))
end
figSave = configureFigureToPDF(figSave);
if eidff(2)
    print(figSave, absBaseFilePath + fileExt(2),'-dpdf','-fillpage')
end
if eidff(3)
    print(figSave, absBaseFilePath + fileExt(3),'-dmeta')
end