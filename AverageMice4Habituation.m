fnOpts = {'UniformOutput', false};
axOpts = {'NextPlot', 'add', 'Color', 'none','Box','off'};
lgOpts = {axOpts{3:end}, 'Location', 'best'};
lnOpts = {'LineStyle', 'none', 'Marker','o','MarkerEdgeColor','none', ...
    'DisplayName'};
figDir = "C:\Users\neuro\seadrive_root\Emilio U\Shared with groups\"+ ...
    "GDrive GrohLab\Projects\00 SC\SC Behaviour\Figures\"+ ...
    "Figure 1\Matlab figures";
% Grouping the movement probability per intensity
miceSelection = [1:3,5];
mmp = arrayfun(@(m) [cat(1, m.Sessions.Intensities), ...
    cat(1, m.Sessions.RollMovProb)], mice(miceSelection), fnOpts{:});
[xu,~,rm] = cellfun(@(x) unique(x(:,1)), mmp, fnOpts{:});
muFlags = cellfun(@(m, u) m(:,1) == u', mmp, xu, fnOpts{:});
% Mean
muVals = cellfun(@(m, f) arrayfun(@(c) mean(m(f(:,c),2)), 1:size(f,2)), ...
    mmp, muFlags, fnOpts{:});
% Quartiles 25, 50, 75
bxPlt = cellfun(@(m, f) arrayfun(@(c) quantile(m(f(:,c),2),3), ...
    1:size(f,2), fnOpts{:}), mmp, muFlags, fnOpts{:});
bxPlt = cellfun(@(x) cat(1, x{:}), bxPlt, fnOpts{:});

x = arrayfun(@(m) arrayfun(@(s) s.Intensities, m.Sessions, fnOpts{:}), ...
    mice(miceSelection), fnOpts{:}); 
while iscell(x); x = cat(1, x{:}); end
y = arrayfun(@(m) arrayfun(@(s) s.RollMovProb, m.Sessions, fnOpts{:}), ...
    mice(miceSelection), fnOpts{:});
while iscell(y); y = cat(1, y{:}); end
%% Average plot
xmu = cat(1, xu{:}); ymu = cat(2, muVals{:})';
[mu_Mdl, muGOF] = fit(xmu, ymu, 'poly1');
miceNames = cat(1, arrayfun(@(x) x.Name, mice(miceSelection), fnOpts{:}));

fig = figure('Name', 'Mean habituation', 'Color', 'w');
ax = axes('Parent', fig, axOpts{:});
scObj = cellfun(@(x,y,n) scatter(ax, x+random(jittNoise, size(x,1), 1), ...
    y', [], 'filled', 'DisplayName', n), xu, muVals, miceNames); 
lObj = plot(mu_Mdl, 'predfun'); set(lObj, 'Color', 'k'); 
lObj(1).DisplayName = 'Trend'; lObj(2).DisplayName = 'Confidence bounds';
lgObj = legend([scObj; lObj(1:2)]); set(lgObj, lgOpts{:})
xlabel(ax, 'Puff intensity [bar]'); ylabel(ax, 'Movement probability')
title(ax, 'Mice movement depending on puff intensity')
saveFigure(fig, fullfile(figDir, "Puff intensity mean dependency"), true, true);
pause;
delete(fig); clearvars fig ax *Obj

%% Quartiles plot
yqr = cellfun(@(x) x(:,2), bxPlt, fnOpts{:}); yqr = cat(1, yqr{:});
[qr_Mdl, qrGOF] = fit(xmu, yqr, 'poly1');

fig = figure('Name', 'Quartiles habituation', 'Color', 'w');
ax = axes('Parent', fig, axOpts{:});
ebObj = cellfun(@(x,y,n,c) errorbar(x+random(jittNoise, size(x,1), 1), ...
    y(:,2), diff(y(:,[1,2]),1,2), diff(y(:,[2,3]),1,2), lnOpts{:} , n, ...
    'MarkerFaceColor', c), xu, bxPlt, miceNames, ...
    mat2cell(clrMap,ones(1,size(clrMap,1)),3)); % converting rows into cells
lObj = plot(qr_Mdl, 'predfun'); set(lObj, 'Color', 'k'); 
lObj(1).DisplayName = 'Trend'; lObj(2).DisplayName = 'Confidence bounds';
lgObj = legend([ebObj; lObj(1:2)]); set(lgObj, lgOpts{:})
xlabel(ax, 'Puff intensity [bar]'); ylabel(ax, 'Movement probability')
title(ax, 'Mice movement depending on puff intensity')
saveFigure(fig, fullfile(figDir, "Puff intensity quartile dependency"), true, true);
pause;
delete(fig); clearvars fig ax *Obj

%% Binning probabilities for all animals
intBinSz = 0.5; intLms = [0, 3]; getRange = @(l, b) l(1):b:l(2);
getBinCent = @(x) mean([x(1:end-1);x(2:end)],1);
intEdge = getRange(intLms, intBinSz); intMemb = discretize(x, intEdge);
intCntr = getBinCent(intEdge);
probHist = arrayfun(@(im) quantile(y(intMemb==im), 3), (1:size(intCntr,2))', ...
    fnOpts{:}); probHist = cat(1, probHist{:}); probHist(isnan(probHist)) = 0;

[bn_Mdl, bnGOF] = fit(intCntr', probHist(:,2), 'poly1');

fig = figure('Name', 'Binning intensities', 'Color', 'w');
ax = axes('Parent', fig, axOpts{:});
hObj = histogram(ax, 'BinEdges',intEdge, 'BinCounts',probHist(:,2), ...
    'EdgeColor','none','FaceColor',0.15*ones(1,3));
erObj = errorbar(ax, intCntr, probHist(:,2), diff(probHist(:,[1,2]),1,2), ...
    diff(probHist(:,[2,3]),1,2), 'LineStyle', 'none', 'Color',0.15*ones(1,3));
lObj = plot(bn_Mdl, 'predfun'); set(lObj, 'Color','k')
legend('off')
xlabel(ax, 'Puff intensity [bar]'); ylabel(ax, 'Movement probability')
title(ax, sprintf('Mice movement depending on puff intensity ( %s)', ...
    sprintf("%d ", miceSelection)))
%%
saveFigure(fig, fullfile(figDir, ...
    sprintf("Puff intensity %.2f binned dependency %smice", intBinSz, ...
    sprintf("%d ", miceSelection))), true, true);
%pause;
delete(fig); clearvars fig ax *Obj