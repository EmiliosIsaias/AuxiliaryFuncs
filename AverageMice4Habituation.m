fnOpts = {'UniformOutput', false};
axOpts = {'NextPlot', 'add', 'Color', 'none','Box','off'};
lgOpts = {axOpts{3:end}, 'Location', 'best'};
lnOpts = {'LineStyle', 'none', 'Marker','o','MarkerEdgeColor','none', ...
    'DisplayName'};
figDir = "C:\Users\neuro\seadrive_root\Emilio U\Shared with groups\"+ ...
    "GDrive GrohLab\Projects\00 SC\SC Behaviour\Figures\"+ ...
    "Figure 1\Matlab figures";
mmp = arrayfun(@(m) [cat(1, m.Sessions.Intensities), ...
    cat(1, m.Sessions.RollMovProb)], mice, fnOpts{:});
[xu,~,rm] = cellfun(@(x) unique(x(:,1)), mmp, fnOpts{:});
muFlags = cellfun(@(m, u) m(:,1) == u', mmp, xu, fnOpts{:});
muVals = cellfun(@(m, f) arrayfun(@(c) mean(m(f(:,c),2)), 1:size(f,2)), ...
    mmp, muFlags, fnOpts{:});
bxPlt = cellfun(@(m, f) arrayfun(@(c) quantile(m(f(:,c),2),3), ...
    1:size(f,2), fnOpts{:}), mmp, muFlags, fnOpts{:});
bxPlt = cellfun(@(x) cat(1, x{:}), bxPlt, fnOpts{:});

%% Average plot
xmu = cat(1, xu{:}); ymu = cat(2, muVals{:})';
[mu_Mdl, muGOF] = fit(xmu, ymu, 'poly1');
miceNames = cat(1, arrayfun(@(x) x.Name, mice, fnOpts{:}));

hf = figure('Name', 'Mean habituation', 'Color', 'w');
ax = axes('Parent', hf, axOpts{:});
scObj = cellfun(@(x,y,n) scatter(ax, x+random(jittNoise, size(x,1), 1), y', [], ...
    'filled', 'DisplayName',n), xu, muVals, miceNames); 
lObj = plot(mu_Mdl, 'predfun'); set(lObj, 'Color', 'k'); 
lObj(1).DisplayName = 'Trend';
lgObj = legend([scObj; lObj(1:2)]); set(lgObj, lgOpts{:})
xlabel(ax, 'Puff intensity [bar]'); ylabel(ax, 'Movement probability')
title(ax, 'Mice movement depending on puff intensity')
saveFigure(hf, fullfile(figDir, "Puff intensity mean dependency"), true, true);
pause;
delete(hf); clearvars hf ax *Obj

%% Quartiles plot
yqr = cellfun(@(x) x(:,2), bxPlt, fnOpts{:}); yqr = cat(1, yqr{:});
[qr_Mdl, qrGOF] = fit(xmu, yqr, 'poly1');

qf = figure('Name', 'Quartiles habituation', 'Color', 'w');
ax = axes('Parent', qf, axOpts{:});
ebObj = cellfun(@(x,y,n,c) errorbar(x+random(jittNoise, size(x,1), 1), ...
    y(:,2), y(:,1), y(:,3), lnOpts{:} , n, 'MarkerFaceColor', c), xu, ...
    bxPlt, miceNames, mat2cell(clrMap,ones(1,size(clrMap,1)),3));
lObj = plot(qr_Mdl, 'predfun'); set(lObj, 'Color', 'k'); 
lObj(1).DisplayName = 'Trend';
lgObj = legend([ebObj; lObj(1:2)]); set(lgObj, lgOpts{:})
xlabel(ax, 'Puff intensity [bar]'); ylabel(ax, 'Movement probability')
title(ax, 'Mice movement depending on puff intensity')
saveFigure(qf, fullfile(figDir, "Puff intensity quartile dependency"), true, true);
pause;
delete(qf); clearvars qf ax *Obj