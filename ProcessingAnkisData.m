% Processing Anki's data
%#ok<*AGROW,*SAGROW> 
% dataPttrn = ...
%     "Z:\Emilio\SuperiorColliculusExperiments\" + ...
%     "Roller\Batch*_beh\WT*\*\*bar*";
% dataPttrn = "Z:\Emilio\SuperiorColliculusExperiments\Roller\"+...
%     "Batch13_beh\WT*\*\*bar*";
dataPttrn = "Z:\Emilio\SuperiorColliculusExperiments\Roller\" + ...
    "Batch13_beh\WT*\2303*\*bar*";
barFlds = dir(dataPttrn);
fn = @(x) fullfile(x.folder, x.name);

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
                struct('Date',currSess,'Intensities',[],'BehIndex',[]);
        else
            mice(mc).Sessions = [mice(mc).Sessions; ...
                struct('Date',currSess,'Intensities',[],'BehIndex',[])];
        end
        sc = sc + 1;
    end
    try
        [outStr, behFigDir] = analyseBehaviour(fn(cf), 'verbose', false,...
            'showPlots', false);
    catch ME
        %dbstop in ProcessingAnkisData.m at 39
        fprintf(1, "Something went wrong with %s\n", fn(cf));
        continue
    end
    consCondNames = arrayfun(@(c) string(c.ConditionName), outStr);
    biFigPttrn = "BehIndex%s";
    biFigPttrn = sprintf(biFigPttrn, sprintf(" %s (%%.3f)", consCondNames));
    [pAreas, ~, behAreaFig] = createBehaviourIndex(outStr);
    outStr = arrayfun(@(bs, ba) setfield(bs,'BehIndex', ba), outStr, pAreas);
    set(behAreaFig, 'UserData', outStr)

    biFN = sprintf(biFigPttrn, pAreas);
    saveFigure(behAreaFig, fullfile(behFigDir, biFN), true);
    mice(mc).Sessions(sc).Intensities = ...
        [mice(mc).Sessions(sc).Intensities; puffInt];
    mice(mc).Sessions(sc).BehIndex = ...
        [mice(mc).Sessions(sc).BehIndex; ...
        outStr.BehIndex];
    close('all')
end

%% Plotting data all together

piFig = figure("Name","Puff intensity vs Movement", "Color","w");
ax = axes("Parent",piFig,"NextPlot","add");
fnOpts = {'UniformOutput', false};
lgOpts = {'Box', 'off', 'Location', 'best', 'NumColumnsMode', 'auto'};
clrMap = lines(size(mice,1)); x = []; y = x; scObj = gobjects(size(mice));
jittNoise = makedist("Normal","mu",0,"sigma",0.025);
for mc = 1:size(mice,1)
    for sc = 1:size(mice(mc).Sessions, 1)
        x = [x; mice(mc).Sessions(sc).Intensities]; 
        y = [y; mice(mc).Sessions(sc).BehIndex];
        scObj(mc) = scatter(ax, ...
            random(jittNoise,size(mice(mc).Sessions(sc).Intensities)) + ...
            mice(mc).Sessions(sc).Intensities, ...
            random(jittNoise,size(mice(mc).Sessions(sc).Intensities))*0+...
            mice(mc).Sessions(sc).BehIndex, ".", ...
            "SizeData", 160, "MarkerEdgeColor", clrMap(mc,:), ...
            "MarkerEdgeAlpha", 0.5);
    end
end
xlabel("Puff intensity [bar]"); ylabel("Movement probability");
set(ax, "Box", "off", "Color", "none");
%[mdl, ~, ~] = polyfit(x, y, 1);
mdl = fit_poly(x, y, 1);
lObj = plot(ax, [0;3], ([0;3].^[1,0])*mdl, "k");
lgObj = legend([scObj; lObj], ...
    cat(1, arrayfun(@(x) x.Name, mice, "UniformOutput",false), 'Trend'));
set(lgObj, "Box", "off", "Location", "best", "NumColumnsMode", "auto")
%% Plotting data. One by one
figDir = "Z:\Emilio\SuperiorColliculusExperiments\Roller\Batch13_beh";
piFig = figure("Name","Puff intensity vs Movement", "Color","w");
clrMap = lines(size(mice,1)); x = []; y = x; scObj = gobjects(size(mice));
jittNoise = makedist("Normal","mu",0,"sigma",0.025);
ax = axes("Parent",piFig,"NextPlot","add");
for mc = 1:size(mice,1)
    for sc = 1:size(mice(mc).Sessions, 1)
        x = [x; mice(mc).Sessions(sc).Intensities]; 
        y = [y; mice(mc).Sessions(sc).BehIndex];
        scObj(mc) = scatter(ax, ...
            random(jittNoise,size(mice(mc).Sessions(sc).Intensities)) + ...
            mice(mc).Sessions(sc).Intensities, ...
            random(jittNoise,size(mice(mc).Sessions(sc).Intensities))*0+...
            mice(mc).Sessions(sc).BehIndex, ".", ...
            "SizeData", 108, "MarkerEdgeColor", clrMap(mc,:));
    end
    xlabel("Puff intensity [bar]"); ylabel("Movement probability");
    lgObj = legend(scObj(mc), mice(mc).Name); set(lgObj, lgOpts{:})
    axis(ax, [-0.5, 3.5, 0, 0.7])
    saveFigure(piFig, fullfile(figDir, sprintf("%s puff vs prob", ...
        mice(mc).Name)), true, true)
    delete(get(ax, "Children"))
end
%%
xlabel("Puff intensity [bar]"); ylabel("Movement probability");
set(ax, "Box", "off", "Color", "none");
[mdl, S, mu] = polyfit(x, y, 1);
lObj = plot(ax, [0;3],([0;3].^[0,1])*mdl', "k");
lgObj = legend([scObj; lObj], ...
    cat(1, arrayfun(@(x) x.Name, mice, "UniformOutput",false), 'Trend'));
set(lgObj, "Box", "off", "Location", "best", "NumColumnsMode", "auto")