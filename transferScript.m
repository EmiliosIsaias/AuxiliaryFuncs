fnOpts = {'UniformOutput', false};
axOpts = {'Box','off','Color','none'};
lgOpts = cat(2, axOpts{1:2}, {'Location','best'});
ffle = @(x) fullfile(x.folder, x.name);
vFiles = dir(fullfile(behDir, "roller*.avi"));
fr = VideoReader(ffle(vFiles));
fr = fr.FrameRate;
rpFiles = dir(fullfile(behDir, "Roller_position*.csv"));
rp = readRollerPositionsFile(ffle(rpFiles));
[vf, rollTx] = getRollerSpeed(rp, fr);
pTms = getCSVTriggers(fullfile(behDir, "Puff.csv")) - rollTx(1);
lTms = getCSVTriggers(fullfile(behDir, "Laser.csv")) - rollTx(1);
trialID = readCSVtimeStamps(fullfile(behDir, "PL.csv"));
trialFlag = false(size(pTms));

condNames = ["Control puff", "Laser100msPuff", "Laser"];

chCond = [1,2];
Nccond = length(chCond);
consCondNames = cellstr(condNames(chCond));
trialFlag(sub2ind(size(pTms), 1:size(pTms,1), ...
    trialID(any(trialID(:,1) == chCond, 2), 1)')) = true;

bvWin = [-0.25, 0.5];
brWin = [1, 400]*1e-3;

[~, vStack] = getStacks(false, round(pTms*fr), 'on', bvWin, fr, fr, ...
    [], vf);
[~, Nt, Na] = size(vStack);
stMdl = fit_poly([1, Nt], bvWin, 1); behTx = ((1:Nt)'.^[1,0]) * stMdl;
bsFlag = behTx < 0;
brFlag = behTx < brWin; brFlag = xor(brFlag(:,1),brFlag(:,2));
% figure; plot(stTx, squeeze(vStack(:,:,trialFlag(:,1))))
%%
figureDir = fullfile(behDir, "BehFigures");
if ~exist(figureDir, "dir")
    if ~mkdir(figureDir)
        error("Could not create figure directory!\n")
    end
end
sSig = squeeze(std(vStack(:,bsFlag,:), [], 2));
[~, sigOrd] = sort(sSig, "descend");
% A bit arbitrary threshold, but enough to remove running trials
sigTh = 2.5; excFlag = false(Na,1);
% excFlag = sSig > sigTh;
ptOpts = {"Color", 0.7*ones(1,3), "LineWidth", 0.2;...
    "Color", "k", "LineWidth",  1.5};
spTh = {0.1:0.1:3}; % Speed threshold
gp = zeros(Nccond, 1, 'single');
rsPttrn = "%s roller speed VW%.2f - %.2f s RM%.2f - %.2f ms EX%d";
pfPttrn = "%s move probability %.2f RW%.2f - %.2f ms EX%d";
rsSgnls = cell(Nccond, 1); mvFlags = cell(Nccond,1); mvpt = mvFlags;
qSgnls = rsSgnls; mat2ptch = @(x) [x(1:end,:)*[1;1]; x(end:-1:1,:)*[1;-1]];
getThreshCross = @(x) sum(x)/size(x,1);
xdf = arrayfun(@(x) ~excFlag & trialFlag(:,x), 1:Nccond, ...
    fnOpts{:});  xdf = cat(2, xdf{:});

for ccond = 1:Nccond
    sIdx = xdf(:,ccond);
    % % Plot speed signals
    fig = figure("Color", "w");
    Nex = sum(xor(sIdx, trialFlag(:,ccond)));
    rsFigName = sprintf(rsPttrn,consCondNames{ccond}, bvWin,...
        brWin*1e3, Nex);
    % Plot all trials
    plot(behTx, squeeze(vStack(:,:,sIdx)), ptOpts{1,:}); hold on;
    % Plot mean of trials
    % Standard deviation
    %rsSgnls{ccond} = [squeeze(mean(vStack(:,:,sIdx),3))',...
    %squeeze(std(vStack(:,:,sIdx),1,3))'];
    % S.E.M.
    rsSgnls{ccond} = [squeeze(mean(vStack(:,:,sIdx),3))',...
        squeeze(std(vStack(:,:,sIdx),1,3))'./sqrt(sum(sIdx))];
    qSgnls{ccond} = quantile(squeeze(vStack(:,:, sIdx)), 3, 2);
    lObj = plot(behTx, rsSgnls{ccond}(:,1), ptOpts{2,:});
    lgnd = legend(lObj,string(consCondNames{ccond}));
    set(lgnd, "Box", "off", "Location", "best")
    set(gca, axOpts{:})
    title(['Roller speed ',consCondNames{ccond}])
    xlabel("Time [s]"); ylabel("Roller speed [cm/s]"); xlim(bvWin)
    saveFigure(fig, fullfile(figureDir, rsFigName), 1)
    % Probability plots
    mvpt{ccond} = getMaxAbsPerTrial(squeeze(vStack(:,:,sIdx)), ...
        brWin, behTx);
    mvFlags{ccond} = compareMaxWithThresh(mvpt{ccond}, spTh);
    gp(ccond) = getAUC(mvFlags{ccond});
    pfName = sprintf(pfPttrn, consCondNames{ccond}, gp(ccond),...
        brWin*1e3, Nex);
    fig = plotThetaProgress(mvFlags(ccond), spTh,...
        string(consCondNames{ccond}));
    xlabel("Roller speed \theta [cm/s]");
    title(sprintf("Trial proportion crossing \\theta: %.3f", gp(ccond)))
    saveFigure(fig, fullfile(figureDir, pfName), 1)
end
clMap = lines(Nccond);
phOpts = {'EdgeColor', 'none', 'FaceAlpha', 0.25, 'FaceColor'};
% Plotting mean speed signals together
fig = figure("Color", "w"); axs = axes("Parent", fig, "NextPlot", "add");
arrayfun(@(x) patch(axs, behTx([1:end, end:-1:1]),...
    mat2ptch(rsSgnls{x}), 1, phOpts{:}, clMap(x,:)), 1:Nccond); hold on
lObj = arrayfun(@(x) plot(axs, behTx, rsSgnls{x}(:,1), "Color", clMap(x,:),...
    "LineWidth", 1.5, "DisplayName", consCondNames{x}), 1:Nccond);
xlabel(axs, "Time [s]"); xlim(axs, bvWin); ylabel(axs, "Roller speed [cm/s]")
set(axs, axOpts{:}); title(axs, "Roller speed for all conditions")
lgnd = legend(axs, lObj); set(lgnd, lgOpts{:})
rsPttrn = "Mean roller speed %s VW%.2f - %.2f s RM%.2f - %.2f ms EX%s SEM";
Nex = sum(trialFlag) - sum(xdf);
rsFigName = sprintf(rsPttrn, sprintf('%s ', consCondNames{:}), bvWin,...
    brWin*1e3, sprintf('%d ', Nex)); saveFigure(fig, fullfile(figureDir, rsFigName), 1)

% Plotting median speed signals together
fig = figure("Color", "w"); axs = axes("Parent", fig, "NextPlot", "add");
arrayfun(@(x) patch(axs, behTx([1:end, end:-1:1]),...
    mat2ptch(qSgnls{x}), 1, phOpts{:}, clMap(x,:)), 1:Nccond); hold on
lObj = arrayfun(@(x) plot(axs, behTx, qSgnls{x}(:,1), "Color", clMap(x,:),...
    "LineWidth", 1.5, "DisplayName", consCondNames{x}), 1:Nccond);
xlabel(axs, "Time [s]"); xlim(axs, bvWin); ylabel(axs, "Roller speed [cm/s]")
set(axs, axOpts{:}); title(axs, "Roller speed for all conditions")
lgnd = legend(axs, lObj); set(lgnd, lgOpts{:})
rsPttrn = "Median roller speed %s VW%.2f - %.2f s RM%.2f - %.2f ms EX%s SEM";
Nex = sum(trialFlag) - sum(xdf);
rsFigName = sprintf(rsPttrn, sprintf('%s ', consCondNames{:}), bvWin,...
    brWin*1e3, sprintf('%d ', Nex)); 
saveFigure(fig, fullfile(figureDir, rsFigName), 1)

% Plotting movement threshold crossings
fig = figure("Color", "w"); axs = axes("Parent", fig, "NextPlot", "add");
mvSgnls = cellfun(getThreshCross, mvFlags, fnOpts{:});
mvSgnls = cat(1, mvSgnls{:}); mvSgnls = mvSgnls';
plot(axs, spTh{1}, mvSgnls);
ccnGP = cellfun(@(x, y) [x, sprintf(' AUC%.3f',y)], consCondNames', ...
    num2cell(gp), fnOpts{:});
lgnd = legend(axs, ccnGP); set(axs, axOpts{:})
set(lgnd, lgOpts{:}); ylim(axs, [0,1])
xlabel(axs, "Roller speed \theta [cm/s]"); ylabel(axs, "Trial proportion")
title(axs, "Trial proportion crossing \theta")
pfPttrn = "Move probability %sRW%.2f - %.2f ms";
pfName = sprintf(pfPttrn, sprintf('%s ', ccnGP{:}), brWin*1e3);
saveFigure(fig, fullfile(figureDir, pfName), 1)
% Tests for movement
prms = nchoosek(1:Nccond,2);
getDistTravel = @(x) squeeze(sum(abs(vStack(:, brFlag, xdf(:,x))), 2));
dstTrav = arrayfun(getDistTravel, 1:Nccond, fnOpts{:});
[p, h, stats] = arrayfun(@(x) ranksum(dstTrav{prms(x,1)}, ...
    dstTrav{prms(x,2)}), 1:size(prms,1), fnOpts{:});
