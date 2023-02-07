%% Auxiliary variables
fnOpts = {'UniformOutput', false};
vars2load = {'gclID', 'popEffects', 'Results', 'logPSTH', ...
    'consCondNames', 'configStructure', 'Counts'};
% Extend file name from file structure
exfn = @(x) fullfile(x.folder, x.name);
getCondNames = @(y) arrayfun(@(x) string(x.ConditionName), y);
removeZeroes = @(c) cellfun(@(y,z) y(z,:), c,...
    cellfun(@(x) ~all(~x,2), c, fnOpts{:}), fnOpts{:});
findStr = @(x, y) contains(x, y, "IgnoreCase", true);
parentDir = "Z:\Emilio\SuperiorColliculusExperiments\"+ ...
    "Roller\Batch11_ephys.MC";
ephFiles = dir(fullfile(parentDir,"*","*","*", ...
    "ephys*\Results\Res*.mat"));
configPrms = cell2mat(arrayfun(@(x) cell2mat(textscan(x.name, ...
    "Res VW%.2f-%.2f ms RW%.2f-%.2f ms SW%.2f-%.2f ms PuffAll.mat")), ...
    ephFiles, fnOpts{:}))*1e-3;
vWins = configPrms(:,1:2); rWins = configPrms(:,3:4);
sWins = configPrms(:,5:6);
% User input, weakest point in the script
mouseTypePoss = ["ChR2"; "eOPN3"]; mtRows = [3;2];
avWins = [-0.3, 0.4;-1.1, 0.4]; arWins = [0.02, 0.05];
asWins = [-0.15, -0.12; -1.05, -1.02];

mouseType = arrayfun(@(x) arrayfun(@(y) findStr(x.folder, y), ...
    mouseTypePoss), ephFiles, fnOpts{:}); mouseType = cat(2, mouseType{:});
prmsFlags = all([ismember(vWins, avWins), ismember(rWins,arWins), ...
    ismember(sWins, asWins)], 2);
Nty = size(mouseTypePoss, 1); Nx = sum(mouseType, 2);
oldCmt = false(Nty,1);
outMat = arrayfun(@(ct) zeros(Nx(ct), mtRows(ct)), 1:Nty, fnOpts{:});
outMat2 = arrayfun(@(ct) zeros(Nx(ct), mtRows(ct)-1), 1:Nty, fnOpts{:});
condNamesCell = cell(Nty,1);
cts = ones(Nty, 1); lc = 0; cnc = 1;
%% Loop through files
for cf = ephFiles'
    lc = lc + 1;
    if ~prmsFlags(lc)
        % Doesn't fulfill the parameter criteria
        continue
    end
    % summary file name
    rfn = exfn(cf); fprintf(1, "Processing %s...\n", rfn)
    mfObj = matfile(rfn); vrsInFile = fieldnames(mfObj);
    % Removing 'Properties' from the variables to load
    if ~all(ismember(vrsInFile(2:end), vars2load))
        fprintf(1, "%s\ndoesn't have all needed variables!\n", rfn)
        continue
    end
    cmt = mouseType(:,lc);
    if ~all(vWins(lc,:) == avWins(cmt,:)) || ...
            ~all(sWins(lc,:) == asWins(cmt,:))
        fprintf(1, "The configuration parameters for this analysis"+...
            " aren't adecuate. Skipping.\n")
        disp(mfObj.configStructure)
        continue
    end
    % Load variables
    cStrct = mfObj.configStructure; condNames = cStrct.ConsideredConditions;
    muMI = mfObj.popEffects; Counts = mfObj.Counts; rwFr = cellfun(...
        @(x) mean(x, "all"), Counts)./ diff(cStrct.Response_window_s);
    ctrlSub = find(findStr(condNames, "Control"));
    lsrSub = find(findStr(condNames, "Delay"));
    fcSub = find(findStr(condNames(lsrSub), "+"));
    if isempty(lsrSub)
        fprintf(1,'No laser condition found in this experiment... ')
        fprintf(1, 'Skipping\n')
        continue
    end
    flpFlag = false;
    if fcSub == 1
        lsrSub = flip(lsrSub); flpFlag = true;
    end
    aux = rwFr([ctrlSub, lsrSub],2);
    outMat{cmt}(cts(cmt),:) = reshape(aux, 1, length(aux));
    cPrm = [repmat(ctrlSub, size(lsrSub,2), 1), lsrSub'];
    [aux2, aux] = ismember(muMI(:,1:2), cPrm, "rows"); 
    MIrows = find(aux2); aux(~aux2) = [];
    if length(aux) == length(lsrSub)
        outMat2{cmt}(cts(cmt),:) = ~flpFlag*muMI(MIrows(aux),3)' + ...
            flpFlag*flip(muMI(MIrows(aux),3))';
    else
        outMat2{cmt}(cts(cmt),:) = ~flpFlag*muMI(1:mtRows(cmt)-1,3)' + ...
            flpFlag*flip(muMI(1:mtRows(cmt)-1,3))';
    end
    if any(cmt~=oldCmt)
        condNamesCell{cnc} = condNames([ctrlSub, lsrSub]);
        oldCmt = cmt; cnc = cnc + 1;
    end
    cts(cmt) = cts(cmt) + 1;
    fprintf(1, '%s mouse, %d laser conditions\n', mouseTypePoss(cmt), ...
        numel(lsrSub))
end
outMat = removeZeroes(outMat); 
outMat = cellfun(@(x) x./x(:,1), outMat, fnOpts{:});
outMat2 = removeZeroes(outMat2);
[Nx, Nc] = cellfun(@size, outMat);
mdPts = cellfun(@median, outMat, fnOpts{:});
%% Plotting results
%figDir = "C:\Users\AcqPc\seadrive_root\Emilio U\Shared with groups\"+ ...
figDir = "C:\Users\neuro\seadrive_root\Emilio U\Shared with groups\"+ ...
    "GDrive GrohLab\Projects\00 SC\SC Behaviour\"+ ...
    "Figures\Figure 2\Matlab figures";
fg2fn = @(x) fullfile(figDir, x);
sctrOpts = {'filled', 'MarkerFaceColor'};
axOpts = {'NextPlot', 'add', 'Color', 'none','Box','off'};
plOpts = {'LineWidth', 0.1, 'Color', 0.3*ones(1,3)};
clrMps = zeros(max(Nc),3,Nty); clrMps(1,:,:) = repmat([0,0,0],1, 1, Nty);
clrFns = {@blues, @greens, @reds};
figs = gobjects(Nty, 1); hfigs = figs;
fgPttrn = "%s ephys MI comparison %d conditions";
[mnLmt, mxLmt] = cellfun(@(x) bounds(x, "all"), outMat2, fnOpts{:});
lmts = [min(cellfun(@min, mnLmt)), max(cellfun(@max, mxLmt))];
miBins = 32; miDelt = diff(lmts)/miBins;
miEdge = linspace(lmts(1)-miDelt, lmts(2)+miDelt, miBins+1);
for ct = 1:Nty
    figs(ct) = figure('Name', mouseTypePoss(ct), 'Color', 'w');
    hfigs(ct) = figure('Name', "Histogram "+mouseTypePoss(ct),'Color','w');
    aux = clrFns{ct}((Nc(ct))*2);
    clrMps(2:Nc(ct),:,ct) = aux(2:Nc(ct),:);
    cax = axes('Parent', figs(ct), axOpts{:});
    plot(cax, outMat{ct}', plOpts{:}); 
    xlim(cax, [0, Nc(ct)+1]);
    permSubs = nchoosek(1:Nc(ct),2);
    for cc = 1:size(outMat{ct},2)
        scatter(cax, cc+zeros(Nx(ct), 1), outMat{ct}(:,cc), ...
            sctrOpts{:}, squeeze(clrMps(cc,:,ct)))
        if cc < size(outMat{ct},2)
            cax2 = subplot(Nc(ct)-1, 1, cc, axOpts{:});
            title(cax2,join(string(condNamesCell{ct}(permSubs(cc,:)))))
            h = histogram(cax2, outMat2{ct}(:,cc), miEdge, ...
                'Normalization','probability', 'EdgeColor', 'none', ...
                'FaceColor',clrMps(cc+1,:,ct)); aux = median(outMat2{ct}(:,cc));
            xline(cax2, aux, '--', 'Color', ...
                0.15*ones(1,3)); 
            text(cax2, aux, max(h.Values), num2str(aux), ...
                "VerticalAlignment", "bottom", "HorizontalAlignment","left")
        end
    end
    xlabel(cax2, "Modulation index"); ylabel(cax2, "Proportion")
    % Plotting the median
    plot(cax, mdPts{ct}, 'k', plOpts{1}, 4); 
    scatter(cax, 1:Nc(ct), mdPts{ct}./mdPts{ct}(1), 72, ...
        clrMps(1:Nc(ct),:,ct), sctrOpts{:}, 'flat');
    xticks(cax, 1:Nc(ct)); xticklabels(cax, condNamesCell{ct,:})
    ylabel(cax, "Relative response");
    title(cax, sprintf("%s relative population response", mouseTypePoss(ct)))
    % Overwriting figures
    saveFigure(figs(ct), fg2fn(sprintf(fgPttrn, mouseTypePoss(ct), Nc(ct))), ...
        true, true)
    saveFigure(hfigs(ct), fg2fn(sprintf("%s modulation index %d conditions", ...
        mouseTypePoss(ct), Nc(ct))), true, true)
end