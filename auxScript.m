%%
fnOpts = {'UniformOutput', false};
lgOpts = {'Box', 'off', 'Location', 'best','AutoUpdate','off'};
cbOpts = {'AxisLocation', 'in', 'Box', 'off', 'Location', 'west',...
    'Color','w','TickDirection','in'};
axOpts = {'Box', 'off', 'Color', 'none','Clipping','off'};
brOpts = {'EdgeColor','none','FaceColor','flat'};
lnOpts = {'LineStyle','--','Color','w'};
repU = @(x) strrep(x,'_','.');
%%
fsStruct1 = getFirstSpikeInfo(relativeSpkTmsStruct1, configStructure1);
fsStruct2 = getFirstSpikeInfo(relativeSpkTmsStruct2, configStructure2);
fsMu = cell2mat(arrayfun(@(x) [fsStruct1(x).FOStats(:,2);fsStruct2(x).FOStats(:,2)], ...
    1:size(fsStruct1,2), fnOpts{:}));
gclID2 = cellfun(@(x) ['5_' x], gclID2, fnOpts{:});
gclID = [gclID1; gclID2];
[pf1, pE1, pS1, fishTest1] = getResponseProbability(relativeSpkTmsStruct1, configStructure1);
[pf2, pE2, pS2, fishTest2] = getResponseProbability(relativeSpkTmsStruct2, configStructure2);
pf = [pf1;pf2]; pE = [pE1;pE2]; pS = [pS1;pS2];
%%

[PSTHpop, psthTx] = getPSTHfromRelSpkStruct(relativeSpkTmsStruct1, configStructure1);
PSTHexp = getPSTHfromRelSpkStruct(relativeSpkTmsStruct2, configStructure2);
PSTH = cat(1, PSTHpop, PSTHexp); 
timeLapse = configStructure1.Viewing_window_s;
responseWindow = configStructure1.Response_window_s;


%%
cs1 = configStructure1; %cs1.BinSize_s = 1e-3;
[PSTHpop, psthTx, Na1] = getPSTHfromRelSpkStruct(relativeSpkTmsStruct1, cs1);
[PSTHexp, ~, Na2] = getPSTHfromRelSpkStruct(relativeSpkTmsStruct2, cs1);
PSTH = cat(1, PSTHpop, PSTHexp); 
timeLapse = configStructure1.Viewing_window_s;
responseWindow = configStructure1.Response_window_s;
spontanousWindow = configStructure1.Spontaneous_window_s;
respFlag = psthTx >= responseWindow; 
respFlag = xor(respFlag(:,1), respFlag(:,2));
sponFlag = psthTx >= spontanousWindow; 
sponFlag = xor(sponFlag(:,1), sponFlag(:,2));
%%
mP = 0;
prpTh = 0.1:0.01:0.9;
for cp = prpTh
    [rclIdx1, mH1, zH1] = getSignificantFlags(Results1, 'ZProp', cp);
    [rclIdx2, mH2, zH2] = getSignificantFlags(Results2, 'ZProp', cp);
    rclIdx = cat(1, rclIdx1, rclIdx2);
    C = sum(all(rclIdx,2)); T = sum(any(rclIdx,2));
    ctxR = sum(rclIdx);
    fprintf(1,"%.2f: M%d, B%d, T%d, C%d, P%.1f%%\n",cp, ctxR, T,C,100*C/T)
    mP = mP + C/T;
end
mP = mP/length(prpTh);
prpTh = 0.55;
[rclIdx1, mH1, zH1] = getSignificantFlags(Results1, 'ZProp', prpTh);
[rclIdx2, mH2, zH2] = getSignificantFlags(Results2, 'ZProp', prpTh);
rclIdx = [rclIdx1; rclIdx2];
C = sum(all(rclIdx,2)); T = sum(any(rclIdx,2));
ctxR = sum(rclIdx); arcIdx = any(rclIdx,2);
fprintf(1,"%.2f: M%d, B%d, T%d, C%d, P%.1f%%\n",prpTh, ctxR, T,C,100*C/T)
%%
% zPSTH = (PSTH(arcIdx,:,:) - mu(arcIdx))./sig(arcIdx);
%% Significance from zPSTH itself -- NOT GOOD RESULTS!
[~, mu, sig] = zscore(PSTH, 1, [2, 3]);
zPSTH = (PSTH - mu)./sig; zPSTH(isinf(zPSTH)|isinf(zPSTH)) = 0;
zM = squeeze(max(zPSTH(:,respFlag,:),[],2));
alph = 0.95;
nDist = makedist('Normal', "mu", 0, "sigma", 1);
signTh = arrayfun(@(z) fminbnd(@(y) ...
    norm(integral(@(x) nDist.pdf(x), -y, y) - z, 2), -3, 3), alph);
rclIdx = zM > signTh; arcIdx = any(rclIdx,2);
%% t-Test for trial crossing proportion
%p0Vec = 1/N:1/N:0.95;
p0Vec = 1/3;
clMp = turbo(length(p0Vec));
thFig = figure; axs = axes('NextPlot','add', 'Parent',thFig,axOpts{:});
rFig = figure; axs2 = axes('NextPlot','add','Parent',rFig,axOpts{:});
set(rFig, "Colormap", turbo(length(p0Vec)))
scatter(axs2, signMat(:,1), signMat(:,2), 50, ...
    "MarkerEdgeAlpha", 0.8, "MarkerEdgeColor", 0.4*ones(1,3))
c1 = 1;
[signMat1, ~, signMatpt1] = zscoreSignificance(Counts1);
signMat1 = cat(2,signMat1{:});
[~, mH1] = getSignificantFlags(Results1);
[signMat2, ~, signMatpt2] = zscoreSignificance(Counts2);
signMat2 = cat(2,signMat2{:});
[~, mH2] = getSignificantFlags(Results2);
pclIdx = [fishTest1.H; fishTest2.H];
mH = [mH1; mH2]; signMat = [signMat1; signMat2];
for p0 = p0Vec
    tRes1 = proportionTest(signMatpt1, p0);
    tRes2 = proportionTest(signMatpt2, p0);
    tRes = [tRes1; tRes2]; 
    % rclIdx = pclIdx & mH;
    rclIdx = logical(tRes); % & mH;
    C = sum(all(rclIdx,2)); T = sum(any(rclIdx,2));
    ctxR = sum(rclIdx); arcIdx = any(rclIdx,2);
    fprintf(1,"%.2f: M%d, B%d, T%d, C%d, P%.1f%%\n",p0, ctxR, T,C,100*C/T)
    %plot(axs, p0, T, "k.", "MarkerEdgeColor",clMp(c1,:))
    %plot(axs, p0, C, "k+", "MarkerEdgeColor",clMp(c1,:))
    plot(axs, p0, T, "k.")
    plot(axs, p0, C, "k+")
    yyaxis(axs, 'right')
    plot(axs, p0, 100*(C/T),"_", "MarkerEdgeColor",0.4*ones(1,3)); 
    yyaxis(axs, 'left')
    scatter(axs2, signMat(any(rclIdx,2),1), signMat(any(rclIdx,2),2), '.',...
        "MarkerEdgeColor", clMp(c1,:)); c1 = c1 + 1;
end
xlabel(axs,'Proportion \theta');
ylabel(axs,'Total (.) | convergent (+) [# units]')
axs.YAxis(2).Label.String = "Convergence (-) [%]";
axs.YAxis(2).Color = 0.4*ones(1,3); axs.YAxis(2).Limits = [0,100];
title(axs, "Responsive units per \theta")
cb = colorbar(axs2, "eastoutside");
cb.Label.String = 'Proportion \theta';
axis(axs2, [0,1,0,1]); xlabel(axs2, "MC responsiveness"); 
ylabel(axs2, "BC responsiveness")
title(axs2, "Responsiveness per \theta")
%% z PSTH, Final
oSel = 2; % Cortex selection 1- MC; 2- BC
[~, ordSubs] = sort(signMat(:,oSel), "descend");
[zPSTH, mu, sig] = zscore(PSTH, 1, [2,3]); % Z-score for all time and both conditions
cLims = [min(zPSTH(:)), max(zPSTH(:))];

zMu = arrayfun(@(x) mean(zPSTH(rclIdx(:,x),:, x)), 1:2 ,fnOpts{:});
zMu = cat   (1, zMu{:});
zLim = round([min(zMu, [], "all"), max(zMu,[], "all")],1);

axsSubs = (12*(0:4))+(1:5)';
figure; axs(1) = subplot(6,12,axsSubs(:),"NextPlot","add"); 
imagesc(axs(1), timeLapse, [], zPSTH(ordSubs(arcIdx(ordSubs)),:,1)); colormap(axs(1), hot);
xlabel(axs(1), "Time [ms]"); title(axs(1), "MC stimulation"); yticks(axs(1), 1:sum(arcIdx));
yticklabels(axs(1), repU(gclID(ordSubs(arcIdx(ordSubs))))); ylabel(axs(1), "Units")
cb = colorbar(axs(1), cbOpts{:}); cb.Label.String = 'Z-score';
plot(axs(1), repmat(responseWindow,2,1), ylim(axs(1)), lnOpts{:})
ylim(axs(1), [1,sum(arcIdx)]+[-1,1]*0.5); xlim(axs(1), timeLapse)
% xticks(axs(1),linspace(timeLapse(1),timeLapse(2),4));
% xticklabels(axs(1), xticks(axs(1))*1e3);
set(get(axs(1),"XAxis"), "Visible", "off")

axsSubs = (12*(0:4))+(7:11)';
axs(2) = subplot(6,12,axsSubs(:), "NextPlot","add");
imagesc(axs(2),timeLapse, [], zPSTH(ordSubs(arcIdx(ordSubs)),:,2)); colormap(axs(2),hot);
xlabel(axs(2), "Time [ms]"); xticklabels(axs(2), xticks(axs(2))*1e3);
title(axs(2), "BC stimulation"); 
plot(axs(2), repmat(responseWindow,2,1), ylim(axs(2)), lnOpts{:})
ylim(axs(2), [1,sum(arcIdx)]+[-1,1]*0.5); xlim(axs(2), timeLapse)
% xticks(axs(2),linspace(timeLapse(1),timeLapse(2),4));
% xticklabels(axs(2), xticks(axs(2))*1e3);
set(get(axs(2),"XAxis"), "Visible", "off")

axsSubs = 6*(1:2:10);
axs(3) = subplot(6,12,axsSubs(:));
barh(axs(3), 1:sum(arcIdx), -rclIdx(ordSubs(arcIdx(ordSubs)),1), "EdgeColor", "none"); 
hold(axs(3), 'on')
barh(axs(3), 1:sum(arcIdx), rclIdx(ordSubs(arcIdx(ordSubs)),2), "EdgeColor", "none");
ylim(axs(3), ylim(axs(1))); %xlabel(axs(3), "Responsive?")
axs(3).XAxis.Visible = 'off';

axsSubs = 12*(1:5);
axs(4) = subplot(6,12,axsSubs); barh(axs(4), 1:sum(arcIdx), signMat(ordSubs(arcIdx(ordSubs)),2))
set(axs(4), axOpts{:}); ylim(axs(4), ylim(axs(1)))
xlabel(axs(4), "Responsivity")
arrayfun(@(x) set(get(x,'YAxis'),'Visible','off'), axs(2:4));
arrayfun(@(x) set(x, "CLim", cLims), axs([1,2]))
axs(4).XAxis.Visible = 'off';

axsSubs = 12*5+(1:5);
axs(5) = subplot(6, 12, axsSubs(:), "NextPlot", "add");
plt = plot(axs(5), psthTx, zMu(1,:)');
xlim(axs(5), timeLapse); 
xticklabels(axs(5), xticks(axs(5))*1e3); xlabel(axs(5), "Time [ms]")
ylabel(axs(5), "Z-score");
ylim(axs(5), zLim)

axsSubs = 12*5+(7:11);
axs(6) = subplot(6,12,axsSubs(:), "NextPlot", "add");
plt = plot(axs(6), psthTx, zMu(2,:)');
xlim(axs(6), timeLapse); 
xticklabels(axs(6), xticks(axs(6))*1e3); xlabel(axs(6), "Time [ms]")
axs(6).YAxis.Visible = 'off';
ylim(axs(6), zLim)

axsSubs = 66;
axs(7) = subplot(6,12,axsSubs);
pie(axs(7), [T-C, C], [0,1])

set(axs, axOpts{:})

%% Z Population PSTH


%%
Nctx = size(PSTH,3);
axs = gobjects(3,1);

ctxName = ["MC", "BC"];
figureDir = "Z:\Berin\0_Paper\Anaesthetised ephys figures";
phPttrn = "%s PSTH %d responsive units VW%.2f-%.2f ms RW%.2f-%.2f ms SW%.2f-%.2f ms (magn %s)";
clMap = [zeros(1,3);ones(1,3)/3];
pclID = gclID(arcIdx);
pfAux = pf(arcIdx,:);
cbLims = [min(zM(arcIdx,:)), max(zM(arcIdx,:))];
oSel = 2; % -------------- CORTEX Selection 1 - MC, 2 - BC
% % % % no specific order
%ordSubs = (1:sum(arcIdx))';
% % % % Order by latency
% [~, ordSubs] = sort(fsMu(arcIdx,oSel));
% % % % Order by magnitude
[~, ordSubs] = sort(zM(arcIdx,oSel));
%clMp = gray(255);
clMp = [0.65*ones(1,3);0,0,0];
mdl = cell2mat(arrayfun(@(x) fit_poly([min(pf(arcIdx,x)),max(pf(arcIdx,x))], ...
    [255, 1], 1), 1:2, fnOpts{:}));
for cctx = 1:size(PSTH,3)
    fig = figure; axs(1) = subplot(6,6,setdiff(1:6*5, 6*(1:5)));
    clSubs = round((pfAux(:,cctx).^[1,0])*mdl(:,cctx));
    psthSig = mean((PSTH(rclIdx(:,cctx),:,cctx)-mu(rclIdx(:,cctx)))./ ...
        sig(rclIdx(:,cctx)), 1, "omitnan");
    imagesc(timeLapse, [], zPSTH(ordSubs,:,cctx)); 
    yticks(axs(1), 1:sum(arcIdx)); colormap(hot)
    yticklabels(axs(1), repU(pclID(ordSubs))); set(axs(1), "CLim", cbLims)
    title(sprintf("%s responsive units", ctxName(cctx)))
    ylabel("Units"); set(get(axs(1), "XAxis"), "Visible", "off")
    cb = colorbar(axs(1), cbOpts{:}); cb.Label.String = "Z-Score";
    cb.Color = ones(1,3); axs(1).NextPlot = 'add';
    line(axs(1), repmat(responseWindow,2,1), ylim(axs(1)), lnOpts{:})
    % mean PSTH
    axs(2) = subplot(6,6,setdiff((6*5+1):6*6, 6^2)); plot(psthTx, psthSig, ...
        "LineWidth", 1.5, "Color", "k"); ylim(axs(2), [-2/3,4]);
    xlim(axs(2),timeLapse); lgnd = legend(axs(2), 'Population PSTH');
    set(lgnd, lgOpts{:}); ylabel(axs(2), 'Z-score')
    xticklabels(axs(2), xticks(axs(2))*1e3); xlabel(axs(2), 'Time [ms]')
    phName = sprintf(phPttrn, ctxName(cctx), ...
        sum(rclIdx(:,cctx)), timeLapse, ...
        responseWindow([1,2,2,1]).*[1,1,-1,-1]*1e3, ctxName(oSel));
    % Bar plot for firing probability
    axs(3) = subplot(6,6,6*(1:5)); b = barh(axs(3),pfAux(ordSubs,cctx), ...
        brOpts{:}); ylim(axs(3), ylim(axs(1)));
    b.CData = clMp(1+(pfAux(ordSubs,cctx)>0.5),:);
    %b.CData = clMp(clSubs(ordSubs),:);
    axs(3).YDir = 'reverse'; axs(3).YAxis.Visible = 'off';
    %yticks(axs(3), 1:sum(arcIdx))
    %yticklabels(axs(3), repU(pclID(ordSubs))); 
    set(axs, axOpts{:});
    linkaxes(axs([1,2]), "x"); linkaxes(axs([1,3]),'y')
    %saveFigure(fig, fullfile(figureDir, phName), 1)
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
line(axs2, repmat(responseWindow,2,1), ylim(axs2), lnOpts{1:3},'k')
%saveFigure(fig2, fullfile(figureDir, phName), 1)
%% 
sponFlag = psthTx >= configStructure1.Spontaneous_window_s;
sponFlag = xor(sponFlag(:,1), sponFlag(:,2));
[~, mu, sig] = zscore(PSTH(:,sponFlag,:),1,2);
zPSTH = (PSTH - mu)./sig; zPSTH(isinf(zPSTH)|isinf(zPSTH)) = 0;
zM = squeeze(max(zPSTH(:,respFlag,:),[],2));
rclIdx = zM > signTh;
arcIdx = any(rclIdx,2);
%%
fig = gobjects(Nctx,1); axs = gobjects(2,1);
for cctx = 1:Nctx
    fig(cctx) = figure("Color", "w"); 
    axs(1) = subplot(1,6,1:5,"Parent", fig(cctx), "NextPlot", "add");
    imagesc(axs(1), timeLapse, [1,sum(arcIdx)], zPSTH(arcIdx,:,cctx)); 
    colormap(hot); hold(axs(1), 'on'); title(axs(1), ctxName(cctx))
    line(axs(1), repmat(responseWindow,2,1), ylim(axs(1)), ...
        "LineStyle","--","Color","w"); xlim(axs(1), timeLapse)
    ylim(axs(1),[1,sum(arcIdx)]+[-1,1]*0.5)
    xticklabels(axs(1), xticks(axs(1))*1e3)
    xlabel(axs(1), "Time [ms]"); ylabel(axs(1), "Units")

    axs(2) = subplot(1,6,6); plot(zM(arcIdx,cctx), 1:sum(arcIdx))
    set(axs(2), "YLim", ylim(axs(1)))
    hold(axs(2),'on');
    line(repmat(signTh,2,1), [1,sum(arcIdx)], "LineStyle", "--", "Color","k")
    set(axs, axOpts{:}); xlabel(axs(2), "Z-score"); 
    ylim(axs(1),[1,sum(arcIdx)]+[-1,1]*0.5)
    axs(2).YAxis.Visible = 'off';
    linkaxes(axs,'y')
end