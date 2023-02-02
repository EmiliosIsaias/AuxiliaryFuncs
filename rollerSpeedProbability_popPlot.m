

getMI = @(x) diff(x, 1, 2)./sum(x, 2);

lvls = 32;
clMap = rdylgn(lvls);
mi2clr = fit_poly([-1, 1], [1, lvls], 1);
phMI = getMI(pharmTable);

muscMean = mean(pharmTable(phMI < 0,:));
ptxMean = mean(pharmTable(phMI > 0, :));
%% Plot
fig = figure("Colormap", clMap);
ax = axes("Parent", fig, "Box", "off", "Color", "none");
miClr = clMap(round(getMI(pharmTable).^[1,0] * mi2clr),:);
%patch(ax, [1,2], [-1, 1])
for cp = 1:size(pharmTable, 1)
    plot(ax, [1,2], pharmTable(cp,:), "Color", miClr(cp,:), ...
        "Marker", ".", "MarkerSize", 15)
    if cp == 1
        hold(ax, "on")
    end
end
pl1 = plot(ax, [1, 2], muscMean, "Color",clMap(1,:), "LineWidth", 3, "DisplayName", ...
    "Muscimol", "Marker",".", "MarkerSize", 20);
pl2 = plot(ax, [1, 2], ptxMean, "Color",clMap(lvls,:), "LineWidth", 3, "DisplayName", ...
    "Picrotoxin", "Marker",".", "MarkerSize", 20);
lgObj = legend([pl1, pl2]); set(lgObj, "Box", "off", "AutoUpdate", "off", ...
    "Location", "best")

boxchart(0.75*ones(size(pharmTable, 1),1), pharmTable(:, 1), ...
    "JitterOutliers", "on", "BoxWidth", 0.25, "Notch", "on", ...
    "BoxFaceColor", 0.5*ones(1,3), "MarkerColor", 0.5*ones(1,3))

boxchart(2.15*ones(sum(phMI<0),1), pharmTable(phMI<0, 2), ...
    "JitterOutliers", "on", "BoxWidth", 0.125, "Notch", "on", ...
    "BoxFaceColor", clMap(1,:), "MarkerColor", clMap(1,:))

boxchart(2.35*ones(sum(phMI>0),1), pharmTable(phMI>0, 2), ...
    "JitterOutliers", "on", "BoxWidth", 0.125, "Notch", "on", ...
    "BoxFaceColor", clMap(lvls,:), "MarkerColor", clMap(lvls,:))

xticks(ax, [1,2]); xticklabels(ax, {'Control', ...
    '{\color[rgb]{0.64, 0, 0.14}Muscimol}/{\color[rgb]{0, 0.4078, 0.2157}PTX}'})
ylabel(ax, "Movement probability")

set(ax, "Box", "off", "Color", "none")
cb = colorbar(ax); cb.Label.String = 'Decrease < -- > Increase';
cb.Ticks = 0:0.25:1; cb.TickLabels = -1:0.5:1; cb.Box = "off";

saveFigure(fig,'Z:\Emilio\SuperiorColliculusExperiments\Roller\Pharma effects',1)
% boxchart(0.75, pharmTable(phMI < 0,:))