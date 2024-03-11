%% Roller speed per trial
rwN = 4;
[trialSubs, chTrSbs] = equiliseTrials(delayFlags, excFlag);
Mxrs = max(cellfun(@(x, y) max(x(y)), mvpt, chTrSbs'));
clrMap = bwg_cm(128);
cSb = find(contains(consCondNames,'Control'));
clMapAux = lines(Nccond);
if cSb ~= 1
    clMap = zeros(Nccond, 3); clMap(cSb,:) = 0.25*ones(1,3);
    clMap(setdiff(1:Nccond,cSb),:) = clMapAux(2:end,:);
else
    clMap = clMapAux;
end
behFig = figure('Color', 'w', 'Name', 'Behaviour');
axs = gobjects(Nccond+1,1);
axSbs = (0:rwN-2)';
for cc = 1:Nccond
    axs(cc) = subplot(rwN, Nccond, (axSbs.^[1,0])*[Nccond;cc], ...
        'Parent', behFig);
    imagesc(axs(cc), behTx*1e3, [], squeeze(vStack(:,:,trialSubs(:,cc))/Mxrs)')
    xlabel(axs(cc), 'Time [ms]')
end
ylabel(axs(1), 'Trials')
axs(Nccond+1) = subplot(rwN, Nccond, 1+Nccond*(rwN - 1):Nccond*rwN);
arrayfun(@(x) patch(axs(end), 1e3*behTx([1:end, end:-1:1]),...
    mat2ptch(rsSgnls{x}), 1, phOpts{:}, clMap(x,:)), 1:Nccond); hold on
lObj = arrayfun(@(x) plot(axs(end), 1e3*behTx, rsSgnls{x}(:,1), ...
    "Color", clMap(x,:), "LineWidth", 1.5, "DisplayName", ...
    consCondNames{x}), 1:Nccond);
xlabel(axs(end), "Time [ms]"); xlim(axs(end), 1e3*bvWin);
ylabel(axs(end), "Roller speed [cm/s]")
set(axs, axOpts{:}); title(axs(end), "Roller speed for all conditions")
lgnd = legend(axs(end), lObj); set(lgnd, lgOpts{:})
arrayfun(@(x) set(x, 'CLim', [-1,1], 'Colormap', clrMap), ...
    axs(1:end-1))
arrayfun(@(x, y) title(x, y{1}), axs(1:end-1), consCondNames')
arrayfun(@(x) set(get(x, 'YAxis'), 'Visible','off'), ...
    axs(setdiff(1:(Nccond + 1), [1, Nccond + 1])));
cb = colorbar(axs(Nccond), axOpts{1:2}, 'Location', 'west', ...
    'TickDirection', 'none', 'Ticks', [-0.9,0.9], 'TickLabels', ...
    {'Backward', 'Forward'}, 'AxisLocation', 'out');
rsPttrn = "Roller speed per trial %sVW%.2f - %.2f ms EX%s%s";
rsptName = sprintf(rsPttrn, sprintf('%s ', consCondNames{:}), bvWin*1e3, ...
    sprintf('%d ', Nex), thrshStr); rsptFile = fullfile(figureDir, rsptName);
saveFigure(behFig, rsptFile, 1)
%%
axSbs = (0:rwN-2)';
clrMap = rocket(128);
fPSTH = PSTH ./ reshape(Na*binSz, 1, 1, Nccond);
zPSTH = zscore(fPSTH, 1, [2,3]); zPopPSTH = squeeze(mean(zPSTH, 1));
Mxe = max(zPopPSTH, [], "all"); Mne = min(zPopPSTH, [], "all");
ephysFig = figure('Color', 'w', 'Name', 'Ephys');
axs = gobjects(Nccond+1,1);
for cc = 1:Nccond
    axs(cc) = subplot(rwN, Nccond, (axSbs.^[1,0])*[Nccond;cc], ...
        'Parent', ephysFig);
    imagesc(axs(cc), timeLapse*1e3, [], zPSTH(:,:,cc), [Mne, Mxe])
    xlabel(axs(cc), 'Time [ms]'); yticks(axs(cc), []);
    title(axs(cc), consCondNames{cc}); colormap(rocket(128))
end
ylabel(axs(1), 'Units');
axs(Nccond+1) = subplot(rwN, Nccond, 1+Nccond*(rwN - 1):Nccond*rwN, ...
    'NextPlot', 'add');
lObj = arrayfun(@(x) plot(axs(end), 1e3*psthTx, zPopPSTH(:,x), ...
    "Color", clMap(x,:), "LineWidth", 1.5, ...
    "DisplayName", consCondNames{x}), 1:Nccond);
lgnd = legend(axs(end), lObj); set(lgnd, lgOpts{:})
xlabel(axs(end), 'Time [ms]'); ylabel(axs(end), "Z-score")
arrayfun(@(x) set(get(x, 'YAxis'), 'Visible','off'), ...
    axs(setdiff(1:(Nccond + 1), [1, Nccond + 1])));
set(axs, axOpts{:});
cb = colorbar(axs(Nccond), axOpts{1:2}, 'Location', 'west', ...
    'TickDirection', 'none', 'AxisLocation', 'in', ...
    'Color', 0.85*ones(1,3)); cb.Label.String = 'Z-score';
ephysPttrn = 'Z-score all-units PSTH %s VW%.2f - %.2f ms Ntrials%s';
ephysName = sprintf(ephysPttrn, sprintf('%s ', consCondNames{:}), ...
    timeLapse*1e3, sprintf(' %d', Na));
ephysFile = fullfile(figureDir, ephysName);
saveFigure(ephysFig, ephysFile, 1); clearvars ephys* axs
%% Log inset for ephys and behaviour
combFig = figure('Color', 'w', 'Name','Behaviour + Ephys');
axs(1) = subplot(2, 1, 1, "Parent", combFig, "NextPlot", "add");
tmBinWidth = diff(10.^logPSTH.Log10BinEdges(:));
flPSTH = squeeze(mean(logPSTH.LogPSTH, 1))./tmBinWidth;
lObj = arrayfun(@(x) semilogx(axs(1), 1e3*logPSTH.TimeAxis, flPSTH(:,x), ...
    "Color", clMap(x,:), "LineWidth", 1.5, ...
    "DisplayName", consCondNames{x}), 1:Nccond);
lgnd = legend(lObj); set(lgnd, lgOpts{:});
title("SC Population activity in logarithmic time scale")
ylabel("Firing rate [Hz]")

axs(2) = subplot(2, 1, 2, "Parent", combFig, "NextPlot", "add");
strt = 50*1e-3;
lObj = arrayfun(@(x) semilogx(axs(2), 1e3*behTx(behTx>strt), ...
    rsSgnls{x}(behTx>strt,1), "Color", clMap(x,:), "LineWidth", 1.5, ...
    "DisplayName", consCondNames{x}), 1:Nccond);
lgnd = legend(lObj); set(lgnd, lgOpts{:});
xlabel(axs(2), "Log_{10} Time [ms]"); ylabel(axs(2), "Roller speed [cm/s]")
set(axs, axOpts{:});
arrayfun(@(x) set(get(x, 'XAxis'), 'Scale', 'log'), axs)
combPttrn = "PSTH + Roller speed%s ephysVW%.2f - %.2f behVW%.2f - %.2f";
combName = sprintf(combPttrn, sprintf(" %s", consCondNames{:}), ...
    responseWindow*1e3, strt*1e3, bvWin(2)*1e3);
combFile = fullfile(figureDir, combName);
saveFigure(combFig, combFile, 1)
%%
figure; clrmap = lines(5);
gObj = gobjects(5, 2);
for cm = 1:size(xmice.DataTable, 3)
    x = (random( normDist, size( xmice.DataTable, [1, 2]) ) + ...
        repmat( 1:3, size( xmice.DataTable, 1 ), 1 ) )';
    y = squeeze( xmice.DataTable(:,[2,1,3],cm) )';
    gObj(cm,:) = line( x, y, ...
        'Marker', 'o', 'MarkerEdgeColor', 'none', ...
        'MarkerFaceColor', clrmap(cm,:), 'Color', clrmap(cm,:), ...
        'DisplayName', xmice.MiceNames(cm));
    hold on
end
x = repmat(1:3, 10, 1);
y = BImat(:,[2,1,3]);
boxchart(x(:),  y(:), "Notch", "on", "BoxFaceColor", 0.15*ones(1,3) )
lgObj = legend(gObj(:,1), "Location", "best", "Box", "off", ...
    "Color", "none", "AutoUpdate", "off");
ylim([0,1]); title("Batch 18"); ylabel("Behaviour index"); xticks(1:3);
xticklabels({'Continuous', 'Control', 'Frequency'})
line([1,2], 1.15*[1,1]*max(BImat(:,1:2), [], "all"), 'Color', 'k', ...
    'Marker', '|')
text(1.5, 1.15*max(BImat(:,1:2), [], "all"), ...
    sprintf( "p=%.3f", signrank( behMI3(:,1) ) ), ...
    "HorizontalAlignment", "center", "VerticalAlignment", "bottom")

line(2:3, 1.15*[1,1]*max(BImat(:,[1,3]), [], "all"), 'Color', 'k', ...
    'Marker', '|')
text(2.5, 1.15*max(BImat(:,[1,3]), [], "all"), ...
    sprintf( "p=%.3f", signrank( behMI3(:,2) ) ), ...
    "HorizontalAlignment", "center", "VerticalAlignment", "bottom")

set(gca, "Box", "off", "Color", "none")