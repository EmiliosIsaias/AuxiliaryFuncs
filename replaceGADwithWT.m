nwbDir = "C:\Users\neuro\seadrive_root\Emilio U\Shared with groups\GDrive GrohLab\Projects\00 SC\SC Anatomy\Source data for heiDATA\NWB\Awake";
nwbFiles = dir(fullfile(nwbDir, "*.nwb"));
nwbFNs = {nwbFiles.name};
expandName = @(x) fullfile(x.folder, x.name);
pathHere = @(x) fullfile(nwbFiles(1).folder, x);

newNwbFNs = replaceBetween(nwbFNs, 1, regexpPattern('[0-9]+'), 'WT');

for cn = nwbFNs
    cn_i = ismember(nwbFNs, cn);
    nwbObj = nwbRead(expandName(nwbFiles(cn_i)));
    nwbObj.identifier = replaceBetween(nwbObj.identifier, 1, regexpPattern('[0-9]+'), 'WT');
    nwbExport(nwbObj, pathHere(newNwbFNs(cn_i))); clearvars("nwbObj")
end