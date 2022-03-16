%%
rclIdx1 = getSignificantFlags(Results1, 'ZProp', 0.75);
rclIdx2 = getSignificantFlags(Results2, 'ZProp', 0.75);
[PSTHpop, psthTx] = getPSTHfromRelSpkStruct(relativeSpkTmsStruct1, configStructure1);
PSTHexp = getPSTHfromRelSpkStruct(relativeSpkTmsStruct2, configStructure2);
PSTH = cat(1, PSTHpop, PSTHexp); rclIdx = cat(1, rclIdx1, rclIdx2);
timeLapse = configStructure1.Viewing_window_s;
responseWindow = configStructure1.Response_window_s;
fnOpts = {'UniformOutput', false};
%%
fsStruct1 = getFirstSpikeInfo(relativeSpkTmsStruct1, configStructure1);
fsStruct2 = getFirstSpikeInfo(relativeSpkTmsStruct2, configStructure2);
fsMu = cell2mat(arrayfun(@(x) [fsStruct1(x).FOStats(:,2);fsStruct2(x).FOStats(:,2)], ...
    1:size(fsStruct1,2), fnOpts{:}));
%%
[~, mu, sig] = zscore(PSTH, 1, [2, 3]);
Nctx = size(PSTH,3);
arcIdx = any(rclIdx,2); axs = gobjects(2,1);
lgOpts = {'Box', 'off', 'Location', 'best'};
cbOpts = {'AxisLocation', 'in', 'Box', 'off', 'Location', 'east'};
axOpts = {'Box', 'off', 'Color', 'none','Clipping','off'};
ctxName = ["MC", "BC"];
figureDir = "Z:\Berin\0_Paper\Anaesthetised ephys figures";
phPttrn = "%s PSTH %d responsive units VW%.2f-%.2f ms RW%.2f-%.2f ms SW%.2f-%.2f ms (latency %s)";
clMap = [zeros(1,3);ones(1,3)/3];
% Order by latency
oSel = 1; % 1 - MC, 2 - BC
[~, ordSubs] = sort(fsMu(arcIdx,oSel));
% [~, ordSubs] = sort(fsMu(arcIdx,2));
% Order by magnitude
% [~, ordSubs] = sort();
% [~, ordSubs] = sort();
pclID = gclID(arcIdx);
for cctx = 1:size(PSTH,3)
    fig = figure; axs(1) = subplot(6,1,1:5);
    zPSTH = (PSTH(arcIdx, :, cctx) - mu(arcIdx))./sig(arcIdx);
    psthSig = mean((PSTH(rclIdx(:,cctx),:,cctx)-mu(rclIdx(:,cctx)))./ ...
        sig(rclIdx(:,cctx)), 1, "omitnan");
    imagesc(timeLapse, [], zPSTH(ordSubs,:)); yticks(axs(1), 1:sum(arcIdx))
    colormap(hot)
    yticklabels(axs(1), repU(pclID(ordSubs))); set(axs(1), "CLim", [0,6])
    title(sprintf("%s responsive units", ctxName(cctx)))
    ylabel("Units"); set(get(axs(1), "XAxis"), "Visible", "off")
    
    axs(2) = subplot(6,1,6); plot(psthTx, psthSig, "LineWidth", 1.5,...
        "Color", "k");  linkaxes(axs, "x"); ylim(axs(2), [-2/3,4])
    xlim(timeLapse); lgnd = legend(axs(2), 'Population PSTH');
    set(lgnd, lgOpts{:}); ylabel(axs(2), 'Z-score')
    xticklabels(axs(2), xticks(axs(2))*1e3); xlabel(axs(2), 'Time [ms]')
    set(axs, axOpts{:}); phName = sprintf(phPttrn, ctxName(cctx), ...
        sum(rclIdx(:,cctx)), timeLapse, ...
        responseWindow([1,2,2,1]).*[1,1,-1,-1]*1e3, ctxName(oSel));
    cb = colorbar(axs(1), cbOpts{:}); cb.Label.String = "Z-Score";
    cb.Color = ones(1,3);
    saveFigure(fig, fullfile(figureDir, phName), 1)
    if cctx == 1
        fig2 = figure; axs2 = axes(fig2, "NextPlot", "add", axOpts{:});
    end
    plot(axs2, psthTx, psthSig, "LineWidth", 1.5,...
        "DisplayName",ctxName(cctx),"Color", clMap(cctx,:))
end
xlim(axs2, timeLapse); xticklabels(axs2, xticks(axs2)*1e3)
lgnd = legend(axs2); set(lgnd, lgOpts{:}); xlabel(axs2, "Time [ms]")
title(axs2, "Population PSTH per cortex"); ylabel(axs2, "Z-score")
phPttrn = "Population PSTH for %s VW%.2f - %.2f ms RW%.2f - %.2f ms SW%.2f - %.2f ms";
phName = sprintf(phPttrn, sprintf('%s ', ctxName), timeLapse,...
    responseWindow([1,2,2,1]).*[1,1,-1,-1]*1e3);
saveFigure(fig2, fullfile(figureDir, phName), 1)
%%