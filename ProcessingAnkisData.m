% Processing Anki's data
%#ok<*AGROW,*SAGROW> 
dataPttrn = ...
    "Z:\Emilio\SuperiorColliculusExperiments\" + ...
    "Roller\Batch*_beh\WT*\*\*bar*";
barFlds = dir(dataPttrn);
fn = @(x) fullfile(x.folder, x.name);

% TODO: Get the mouse, date and intensity ordered in a structure array.
oldSess = ""; oldMouse = "";
mice = [];
% Mice (mc) and session (sc) counters
mc = 0; sc = 0;
for cf = barFlds'
    fprintf("Processing %s\n", fn(cf))
    puffInt = str2double(extractBefore(cf.name, "bar"));
    [restPath, currSess] = fileparts(cf.folder);
    [restPath, currMouse] = fileparts(restPath);
    if string(oldMouse) ~= string(currMouse)
        oldMouse = currMouse;
        mice = [mice; struct('Name', currMouse, 'Sessions',[])]; 
        mc = mc + 1;
        sc = 0; oldSess = "";
    end
    if string(oldSess) ~= string(currSess)
        oldSess = currSess;
        if ~isfield(mice, 'Sessions')
            mice(mc).Sessions =...
                struct('Date',currSess,'Intensities',[],'RollMovProb',[]);
        else
            mice(mc).Sessions = [mice(mc).Sessions; ...
                struct('Date',currSess,'Intensities',[],'RollMovProb',[])];
        end
        sc = sc + 1;
    end
    try
        outStr = analyseBehaviour(fn(cf), 'verbose', false,...
            'showPlots', false);
    catch ME
        dbstop in ProcessingAnkisData.m at 39
        fprintf(1, "Som ething went wrong with %s\n", fn(cf));
        continue
    end
    mice(mc).Sessions(sc).Intensities = ...
        [mice(mc).Sessions(sc).Intensities; puffInt];
    mice(mc).Sessions(sc).RollMovProb = ...
        [mice(mc).Sessions(sc).RollMovProb; outStr(4).MovProb];
    close all
end