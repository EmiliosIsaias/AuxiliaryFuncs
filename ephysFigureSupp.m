%%
fnOpts = {'UniformOutput', false};
lgOpts = {'Box', 'off', 'Location', 'best','AutoUpdate','off'};
cbOpts = {'AxisLocation', 'in', 'Box', 'off', 'Location', 'west',...
    'Color','w','TickDirection','in'};
axOpts = {'Box', 'off', 'Color', 'none','Clipping','off'};
brOpts = {'EdgeColor','none','FaceColor','flat'};
lnOpts = {'LineStyle','--','Color','w'};
dgOpts = {'EdgeColor', 'none'};
repU = @(x) strrep(x,'_','.'); m = 1e-3; k = 1/m;
gclID = repU(gclID);
%% Refractory period violation cleaning
% Overwriting the spike structure!
relativeSpkTmsStruct_laser = removeRefPViolationSpkFromStruct(relativeSpkTmsStruct_laser);
relativeSpkTmsStruct_whisker = removeRefPViolationSpkFromStruct(relativeSpkTmsStruct_whisker);

%% First spike in response window and its median and mean
fsStruct_laser = getFirstSpikeInfo(relativeSpkTmsStruct_laser, configStructure_laser);
fsStruct_whisker = getFirstSpikeInfo(relativeSpkTmsStruct_whisker, configStructure_whisker);
fsMu = [cell2mat(arrayfun(@(x) x.MedMeanStd(:,2), fsStruct_laser, ...
    fnOpts{:})) fsStruct_whisker.MedMeanStd(:,2)];
fsMed = [cell2mat(arrayfun(@(x) x.MedMeanStd(:,1), fsStruct_laser, ...
    fnOpts{:})) fsStruct_whisker.MedMeanStd(:,1)];
%% 
[pf_l, pE_l, pS_l, fishTest_l] = ...
    getResponseProbability(relativeSpkTmsStruct_laser, configStructure_laser);
[pf_w, pE_w, pS_w, fishTest_w] = ...
    getResponseProbability(relativeSpkTmsStruct_whisker, configStructure_whisker);
pf = [pf_l, pf_w]; pE = [pE_l, pE_w]; pS = [pS_l, pS_w];

%%
[PSTH_laser, psthTx, Na_l] = getPSTHfromRelSpkStruct(relativeSpkTmsStruct_laser, configStructure_laser);
[PSTH_whisker, ~, Na_w] = getPSTHfromRelSpkStruct(relativeSpkTmsStruct_whisker, configStructure_whisker);
PSTH = cat(3, PSTH_laser, PSTH_whisker); 
timeLapse = configStructure_laser.Viewing_window_s;
responseWindow_laser = configStructure_laser.Response_window_s;
spontanousWindow_laser = configStructure_laser.Spontaneous_window_s;
responseWindow_whisker = configStructure_whisker.Response_window_s;
spontanousWindow_whisker = configStructure_whisker.Spontaneous_window_s;
respFlag_laser = psthTx >= responseWindow_laser; 
respFlag_laser = xor(respFlag_laser(:,1), respFlag_laser(:,2));
respFlag_whisker = psthTx >= responseWindow_whisker; 
respFlag_whisker = xor(respFlag_whisker(:,1), respFlag_whisker(:,2));
% sponFlag = psthTx >= spontanousWindow_laser; 
% sponFlag = xor(sponFlag(:,1), sponFlag(:,2));
%% ISI violations measurement
rpTh = [0.5, 1, 1.5, 2, 5, 10]*m;
arpTh = 1.5*m;
[bc, be] = prepareLogBinEdges([1e-5, 1e2], 128);

pSpkTms = cellfun(@(x) x./fs, pSpkSubs, fnOpts{:});
pIsi = cellfun(@diff, pSpkTms, fnOpts{:});
spurSpkFlags = cellfun(@(x) x<m, pIsi, fnOpts{:});
pSpkTmsClean = cellfun(@(x, y) x([false;~y]), pSpkTms, spurSpkFlags, fnOpts{:});
pIsiClean = cellfun(@diff, pSpkTmsClean, fnOpts{:});
% isiProp = cellfun(@(x) sum(x<arpTh)/(numel(x)+1), pIsi);
% lisiDist = cellfun(@(x) histcounts(x, 10.^be), pIsi, fnOpts{:});
% lisiDist = cellfun(@(x) x./diff(10.^be), lisiDist, fnOpts{:});
% lisiDist = cellfun(@(x) x/sum(x), lisiDist, fnOpts{:});
% arpV = cell2mat(cellfun(@(x) interp1(10.^bc, cumsum(x), rpTh), ...
%    lisiDist, fnOpts{:}));
% figure; imagesc(arpV*1e2, [0, 3]); xticklabels(rpTh*k)
isiProp = cellfun(@(x) sum(x<arpTh)/(numel(x)+1), pIsiClean);
lisiDist = cellfun(@(x) histcounts(x, 10.^be), pIsiClean, fnOpts{:});
lisiDist = cellfun(@(x) x/sum(x), lisiDist, fnOpts{:});
arpV = cell2mat(cellfun(@(x) interp1(10.^bc, cumsum(x), rpTh), ...
    lisiDist, fnOpts{:}));



%% Testing thresholds
mP = 0;
prpTh = 0.1:0.01:0.9;
for cp = prpTh
    [rclIdx_laser, mH_l, zH_l] = getSignificantFlags(Results_laser, 'ZProp', cp);
    [rclIdx_whisker, mH_w, zH_w] = getSignificantFlags(Results_whisker, 'ZProp', cp);
    rclIdx = cat(2, rclIdx_laser, rclIdx_whisker);
    Cm = sum(all(rclIdx(:,[1,3]),2));
    Cb = sum(all(rclIdx(:,[2,3]),2));
    Call = sum(all(rclIdx, 2));
    Tl = sum(any(rclIdx(:,1:2),2));
    Tw = sum(rclIdx_whisker);
    Tw = sum(any(rclIdx,2));
    ctxR = sum(rclIdx);
    fprintf(1,"%.2f: MW%d, BW%d, Call%d, C%d, P%.1f%%\n",cp, ctxR, Tl,Cm,100*Cm/Tl)
    mP = mP + Cm/Tl;
end
mP = mP/length(prpTh);
prpTh = 0.55;
[rclIdx_laser, mH_l, zH_l] = getSignificantFlags(Results_laser, 'ZProp', prpTh);
[rclIdx_whisker, mH_w, zH_w] = getSignificantFlags(Results_whisker, 'ZProp', prpTh);
rclIdx = [rclIdx_laser; rclIdx_whisker];
Cm = sum(all(rclIdx,2)); Tl = sum(any(rclIdx,2));
ctxR = sum(rclIdx); arcIdx = any(rclIdx,2);
fprintf(1,"%.2f: M%d, B%d, T%d, C%d, P%.1f%%\n",prpTh, ctxR, Tl,Cm,100*Cm/Tl)
%% Significance from zPSTH itself -- NOT GOOD RESULTS!
%{
[~, mu, sig] = zscore(PSTH, 1, [2, 3]);
zPSTH = (PSTH - mu)./sig; zPSTH(isinf(zPSTH)|isinf(zPSTH)) = 0;
zM = squeeze(max(zPSTH(:,respFlag_laser,:),[],2));
alph = 0.95;
nDist = makedist('Normal', "mu", 0, "sigma", 1);
signTh = arrayfun(@(z) fminbnd(@(y) ...
    norm(integral(@(x) nDist.pdf(x), -y, y) - z, 2), -3, 3), alph);
rclIdx = zM > signTh; arcIdx = any(rclIdx,2);
%}
%% t-Test for trial crossing proportion
N = 30;
p0Vec = 1/N:1/N:1-1/N;
% p0Vec = 1/3;
clMp = turbo(length(p0Vec));
thFig = figure; axs = axes('NextPlot','add', 'Parent',thFig,axOpts{:});
rFig = figure; axs2 = axes('NextPlot','add','Parent',rFig,axOpts{:});
set(rFig, "Colormap", turbo(length(p0Vec)))

[signMat_laser, ~, signMatpt_l] = zscoreSignificance(Counts_laser);
signMat_laser = cat(2, signMat_laser{:});
[~, mH_l] = getSignificantFlags(Results_laser);
[signMat_whisker, ~, signMatpt_w] = zscoreSignificance(Counts_whisker);
signMat_whisker = cat(2,signMat_whisker{:});
[~, mH_w] = getSignificantFlags(Results_whisker);

mH = [mH_l, mH_w]; signMat = [signMat_laser, signMat_whisker];
scatter3(axs2, signMat(:,1), signMat(:,2), 50, ...
    "MarkerEdgeAlpha", 0.8, "MarkerEdgeColor", 0.4*ones(1,3))
c1 = 1;
tRes_l = proportionTest(signMatpt_l, 1/3);
Tl = sum(any(tRes_l,2));
for p0 = p0Vec
    tRes_w = proportionTest(signMatpt_w, p0);
    tRes = [tRes_l, tRes_w];
    pclIdx = logical(tRes);
    Call = sum(all(pclIdx,2));
    Cm = sum(all(pclIdx(:,[1, 3]), 2));
    Cb = sum(all(pclIdx(:,[2, 3]), 2));
    arcIdx = any(pclIdx,2); Tw = sum(pclIdx(:,3));
    fprintf(1,"%.2f: W:%d, Mc%d (%.2f%%), Bc%d (%.2f%%), Xc%d (%.2f%%)\n",p0, Tw, Cm, Cm/Tl, Cb, Cb/Tl, Call, Call/Tl)
    plot(axs, p0, Tw, "k.")
    plot(axs, p0, Call, "k+")
    yyaxis(axs, 'right')
    plot(axs, p0, 100*(Cm/Tl),"_", "MarkerEdgeColor",0.4*ones(1,3)); 
    yyaxis(axs, 'left')
    scatter3(axs2, signMat(arcIdx,1), signMat(arcIdx,2), signMat(arcIdx,3), '.',...
        "MarkerEdgeColor", clMp(c1,:)); c1 = c1 + 1;
    if Tw >= 1
        figure('Name', sprintf('%.2f, %d', p0, Tw));
        imagesc(timeLapse, [], zPSTH(tRes_w==1,:,3)); colormap(rocket);
        yticks(1:Tw); yticklabels(gclID(tRes_w==1))
    end
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

zMu = arrayfun(@(x) mean(zPSTH(rclIdx(:,x),:, x)), 1:size(zPSTH,3),fnOpts{:});
zMu = cat(1, zMu{:});
zLim = round([min(zMu, [], "all"), max(zMu,[], "all")],1);

axsSubs = (12*(0:4))+(1:5)';
figure; axs(1) = subplot(6,12,axsSubs(:),"NextPlot","add"); 
imagesc(axs(1), timeLapse, [], zPSTH(ordSubs(arcIdx(ordSubs)),:,1)); colormap(axs(1), rocket);
xlabel(axs(1), "Time [ms]"); title(axs(1), "MC stimulation"); yticks(axs(1), 1:sum(arcIdx));
yticklabels(axs(1), repU(gclID(ordSubs(arcIdx(ordSubs))))); ylabel(axs(1), "Units")
cb = colorbar(axs(1), cbOpts{:}); cb.Label.String = 'Z-score';
plot(axs(1), repmat(responseWindow_laser,2,1), ylim(axs(1)), lnOpts{:})
ylim(axs(1), [1,sum(arcIdx)]+[-1,1]*0.5); xlim(axs(1), timeLapse)
% xticks(axs(1),linspace(timeLapse(1),timeLapse(2),4));
% xticklabels(axs(1), xticks(axs(1))*1e3);
set(get(axs(1),"XAxis"), "Visible", "off")

axsSubs = (12*(0:4))+(7:11)';
axs(2) = subplot(6,12,axsSubs(:), "NextPlot","add");
imagesc(axs(2),timeLapse, [], zPSTH(ordSubs(arcIdx(ordSubs)),:,2)); colormap(axs(2),rocket);
xlabel(axs(2), "Time [ms]"); xticklabels(axs(2), xticks(axs(2))*1e3);
title(axs(2), "BC stimulation"); 
plot(axs(2), repmat(responseWindow_laser,2,1), ylim(axs(2)), lnOpts{:})
ylim(axs(2), [1,sum(arcIdx)]+[-1,1]*0.5); xlim(axs(2), timeLapse)
% xticks(axs(2),linspace(timeLapse(1),timeLapse(2),4));
% xticklabels(axs(2), xticks(axs(2))*1e3);
set(get(axs(2),"XAxis"), "Visible", "off")

axsSubs = 6*(1:2:10);
axs(3) = subplot(6,12,axsSubs(:));
barh(axs(3), 1:sum(arcIdx), -rclIdx(ordSubs(arcIdx(ordSubs)),1), dgOpts{:}); 
hold(axs(3), 'on')
barh(axs(3), 1:sum(arcIdx), rclIdx(ordSubs(arcIdx(ordSubs)),2), dgOpts{:});
ylim(axs(3), ylim(axs(1))); %xlabel(axs(3), "Responsive?")
axs(3).XAxis.Visible = 'off';

axsSubs = 12*(1:5);
axs(4) = subplot(6, 12, axsSubs); 
% barh(axs(4), 1:sum(arcIdx), ...
%     signMat(ordSubs(arcIdx(ordSubs)),2), dgOpts{:})
bV = barh(axs(4), isiProp(ordSubs(arcIdx(ordSubs)))*1e2, dgOpts{:});
%bV = barh(axs(4), arpV(ordSubs(arcIdx(ordSubs)), 1:4)*1e2, dgOpts{:});
xl = xline(axs(4), 3, 'LineStyle', '--', 'DisplayName', '3%');
set(axs(4), axOpts{:}); ylim(axs(4), ylim(axs(1)))
% xlabel(axs(4), "Responsivity")
% lgnd = legend(bV, string(k*rpTh)+" ms"); set(lgnd, lgOpts{:});
lgnd = legend(bV, string(k*arpTh)+" ms"); set(lgnd, lgOpts{:});
xlabel(axs(4), "ISI violation [%]")
arrayfun(@(x) set(get(x,'YAxis'),'Visible','off'), axs(2:4));
arrayfun(@(x) set(x, "CLim", cLims), axs([1,2]))
title(axs(4), sprintf("ISI < %.2f ms", arpTh*k))
% text(axs(4), 3, find(isiProp(ordSubs(arcIdx(ordSubs)))*1e2<3,1,'first'), ...
%     '3%')
text(axs(4), 3, find(arpV(ordSubs(arcIdx(ordSubs)),4)*1e2<3,1,'first'), ...
    '3%')
axs(4).XAxis.Visible = 'on';
axs(4).YAxis.Visible = 'off';

axsSubs = 12*5+(1:5);
axs(5) = subplot(6, 12, axsSubs(:), "NextPlot", "add");
plot(axs(5), psthTx, zMu(1,:)');
xlim(axs(5), timeLapse); 
xticklabels(axs(5), xticks(axs(5))*1e3); xlabel(axs(5), "Time [ms]")
ylabel(axs(5), "Z-score");
ylim(axs(5), zLim)

axsSubs = 12*5+(7:11);
axs(6) = subplot(6,12,axsSubs(:), "NextPlot", "add");
plot(axs(6), psthTx, zMu(2,:)'); xlim(axs(6), timeLapse); 
xticklabels(axs(6), xticks(axs(6))*1e3); xlabel(axs(6), "Time [ms]")
axs(6).YAxis.Visible = 'off';
ylim(axs(6), zLim)

axsSubs = 66;
axs(7) = subplot(6,12,axsSubs);
pie(axs(7), [Tl-Cm, Cm], [0,1])

set(axs, axOpts{:})

%% Whisker responding units
oSel = 3;
[~, ordSubs] = sort(fsMed(:,oSel), "descend");
cLims = [min(zPSTH(mH_w, :, 3), [], "all"), ...
    max(zPSTH(mH_w, :, 3), [], "all")];

figure('Name','Whisker responding units', 'Color','w'); 
imagesc(timeLapse*k, [], zPSTH(ordSubs(mH_w(ordSubs)), :, 3));
title('Whisker responsive units'); set(gca, axOpts{:});
colormap(rocket); yticks(1:sum(mH_w)); 
yticklabels(gclID(ordSubs(mH_w(ordSubs)))); ylabel('Units')
xlabel('Time [ms]'); %%xline(5:5:15, 'w--')



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
    line(axs(1), repmat(responseWindow_laser,2,1), ylim(axs(1)), lnOpts{:})
    % mean PSTH
    axs(2) = subplot(6,6,setdiff((6*5+1):6*6, 6^2)); plot(psthTx, psthSig, ...
        "LineWidth", 1.5, "Color", "k"); ylim(axs(2), [-2/3,4]);
    xlim(axs(2),timeLapse); lgnd = legend(axs(2), 'Population PSTH');
    set(lgnd, lgOpts{:}); ylabel(axs(2), 'Z-score')
    xticklabels(axs(2), xticks(axs(2))*1e3); xlabel(axs(2), 'Time [ms]')
    phName = sprintf(phPttrn, ctxName(cctx), ...
        sum(rclIdx(:,cctx)), timeLapse, ...
        responseWindow_laser([1,2,2,1]).*[1,1,-1,-1]*1e3, ctxName(oSel));
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
    responseWindow_laser([1,2,2,1]).*[1,1,-1,-1]*1e3);
line(axs2, repmat(responseWindow_laser,2,1), ylim(axs2), lnOpts{1:3},'k')
%saveFigure(fig2, fullfile(figureDir, phName), 1)
%%
fig = gobjects(Nctx,1); axs = gobjects(2,1);
for cctx = 1:Nctx
    fig(cctx) = figure("Color", "w"); 
    axs(1) = subplot(1,6,1:5,"Parent", fig(cctx), "NextPlot", "add");
    imagesc(axs(1), timeLapse, [1,sum(arcIdx)], zPSTH(arcIdx,:,cctx)); 
    colormap(hot); hold(axs(1), 'on'); title(axs(1), ctxName(cctx))
    line(axs(1), repmat(responseWindow_laser,2,1), ylim(axs(1)), ...
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
%% Firing rates
% Spontaneous
rsFr = arrayfun(@(x) quantile(spFr(rclIdx(:,x)), 3), 1:2, fnOpts{:});
% Evoked
reFr = arrayfun(@(x) squeeze(quantile(mean(PSTHfr(rclIdx(:,x), respFlag_laser, ...
    x), 2),3)), 1:2, fnOpts{:});
% Convergence evoked
ceFr = arrayfun(@(x) quantile(mean(PSTHfr(all(rclIdx,2), respFlag_laser, x), 2), ...
    3), 1:2, fnOpts{:});
% Convergence spontaneous
csFr = quantile(spFr(all(rclIdx,2)), 3);
%% Onset latency
fspk = cat(3,fsStruct_laser.MedMeanStd); fspk = cat(1, fspk, cat(3, fsStruct_whisker.MedMeanStd));
onLat = arrayfun(@(x) quantile(fspk(rclIdx(:,x), 1, x),3, 1)*1e3, 1:2, ...
    fnOpts{:});