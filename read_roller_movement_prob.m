fnOpts = {'UniformOutput', false};
% Extend file name from file structure
exfn = @(x) fullfile(x.folder, x.name);
getCondNames = @(y) arrayfun(@(x) string(x.ConditionName), y);
findStr = @(x, y) contains(x, y, "IgnoreCase", true);
parentDir = "Z:\Emilio\SuperiorColliculusExperiments\Roller\Batch11_ephys.MC";
summFiles = dir(fullfile(parentDir, "*","*","*","Behaviour/Simple summary.mat")); 
mouseTypePoss = ["ChR2"; "eOPN3"]; mtRows = [3;2];

mouseType = arrayfun(@(x) arrayfun(@(y) findStr(x.folder, y), ...
    mouseTypePoss), summFiles, 'UniformOutput', false);
mouseType = cat(2, mouseType{:}); oldCmt = false(Nty,1);

Nty = size(mouseTypePoss, 1); Nx = sum(mouseType, 2);
outMat = arrayfun(@(ct) zeros(Nx(ct), mtRows(ct)), 1:Nty, fnOpts{:});
condNamesCell = cell(Nty,1);
cts = ones(Nty, 1); lc = 0; cnc = 1;
for cf = summFiles'
    % summary file name
    sfn = exfn(cf); fprintf(1, "Processing %s...\n", sfn)
    lc = lc + 1;
    load(sfn, 'summStruct'); condNames = getCondNames(summStruct);
    ctrlSub = find(findStr(condNames, "Control"));
    lsrSub = find(findStr(condNames, "Delay"));
    fcSub = find(findStr(condNames(lsrSub), "+"));
    if isempty(lsrSub)
        fprintf(1,'No laser condition found in this experiment... ')
        fprintf(1, 'Skipping\n')
        continue
    end
    if fcSub == 1
        lsrSub = flip(lsrSub);
    end
    cmt = mouseType(:,lc); aux = arrayfun(@(x) x.Results(4).MovProbability, ...
        summStruct([ctrlSub,lsrSub])); 
    outMat{cmt}(cts(cmt),:) = aux;
    if any(cmt~=oldCmt)
        condNamesCell{cnc} = condNames([ctrlSub, lsrSub]);
        oldCmt = cmt; cnc = cnc + 1;
    end
    cts(cmt) = cts(cmt) + 1; 
    fprintf(1, '%s mouse, %d laser conditions\n', mouseTypePoss(cmt), ...
        numel(lsrSub))
end
outMat = cellfun(@(y,z) y(z,:), outMat, ...
    cellfun(@(x) ~all(~x,2), outMat, fnOpts{:}), fnOpts{:});
[Nx, Nc] = cellfun(@size, outMat);
mdPts = cellfun(@median, outMat, fnOpts{:});
%% Plotting results
figDir = "C:\Users\AcqPc\seadrive_root\Emilio U\Shared with groups\"+ ...
    "GDrive GrohLab\Projects\00 SC\SC Behaviour\"+ ...
    "Figures\Figure 2\Matlab figures";
fg2fn = @(x) fullfile(figDir, x);
sctrOpts = {'filled', 'MarkerFaceColor'};
axOpts = {'NextPlot', 'add', 'Color', 'none','Box','off'};
plOpts = {'LineWidth', 0.1, 'Color', 0.3*ones(1,3)};
clrMps = zeros(max(Nc),3,Nty); clrMps(1,:,:) = repmat([0,0,0],1, 1, Nty);
clrFns = {@blues, @greens, @reds};
figs = gobjects(Nty, 1); 
fgPttrn = "%s movement probability comparison %d conditions";
for ct = 1:Nty
    figs(ct) = figure('Name', mouseTypePoss(ct), 'Color', 'w');
    aux = clrFns{ct}((Nc(ct))*2);
    clrMps(2:Nc(ct),:,ct) = aux(2:Nc(ct),:);
    cax = axes('Parent', figs(ct), axOpts{:});
    plot(cax, outMat{ct}', plOpts{:}); xlim(cax, [0, Nc(ct)+1]);
    for cc = 1:size(outMat{ct},2)
        scatter(cax, cc+zeros(Nx(ct), 1), outMat{ct}(:,cc), ...
            sctrOpts{:}, squeeze(clrMps(cc,:,ct)))
    end
    % Plotting the median
    plot(cax, mdPts{ct}, 'k', plOpts{1}, 4); scatter(cax, 1:Nc(ct), ...
        mdPts{ct}, 72, clrMps(1:Nc(ct),:,ct), sctrOpts{:}, 'flat'); 
    xticks(cax, 1:Nc(ct)); xticklabels(cax, condNamesCell{ct,:})
    ylabel(cax, "Movement probability"); 
    title(cax, sprintf("%s movement probability", mouseTypePoss(ct)))
    % Overwriting figures
    saveFigure(figs(ct), fg2fn(sprintf(fgPttrn, mouseTypePoss(ct), Nc(ct))), ...
        true, true)
end