% Processing Anki's data
dataPttrn = ...
    "Z:\Emilio\SuperiorColliculusExperiments\" + ...
    "Roller\Batch*_beh\WT*\*\*bar*";
barFlds = dir(dataPttrn);
fn = @(x) fullfile(x.folder, x.name);

for cf = barFlds'
    fprintf("Processing %s\n", fn(cf))
    puffInt = str2double(extractBefore(cf.name, "bar"));
    % TODO: Get the mouse, date and intensity ordered in a structure array.
end