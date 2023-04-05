%%
axOpts = {'Box', 'off', 'Color', 'none'};
bxOpts = {'Notch', 'on', 'MarkerStyle', 'none', 'BoxWidth', 0.25, ...
    'BoxFaceColor'}; box_sep = 0.15;
myRng = @(x) range(x, "all");
mvRng = cellfun(myRng, behStack);
fgOpts = {'Color', 'w', 'NextPlot', 'add'};
jDist = makedist('Normal', 'mu', 0, 'sigma', 1/32);
% Normality assumption
nDist = makedist('Normal', "mu", 0, "sigma", 1);
alph = 1 - [0.2, 0.1 , 5e-2, 1e-3];
signTh = arrayfun(@(z) fminbnd(@(y) ...
    norm(integral(@(x) pdf(nDist, x), -y, y) - z, 2), -6, 6), alph);
zeroDiv = @(dv) 1./ (dv.*(dv~=0) + 1.*(dv==0));
myzs = @(x, mu, st, c) (c.*(x - mu)) .* zeroDiv(st);
isWhiskOutlier = @(x, Qs, Qr) (Qs(1,:) - 1.5*Qr) > x(:)' |...
    (Qs(2,:) + 1.5*Qr) < x(:)';
getMI = @(a, b) (b - a) .* zeroDiv(a + b);

sponMu = cellfun(@(bs) mean(bs(bsFlag,:), 1, "omitnan"), behStack, fnOpts{:});
sponSig = cellfun(@(bs) std(bs(bsFlag,:), 0, 1, "omitnan"), behStack, fnOpts{:});

respSig = cellfun(@(bs) std(bs(brFlag,:), 0, 1, "omitnan"), behStack, fnOpts{:});

sponQs = cellfun(@(bs) quantile(bs(bsFlag,:), [1,3]/4), behStack, fnOpts{:});
sponMed = cellfun(@(bs) median(bs(bsFlag,:), 1, "omitnan"), behStack, fnOpts{:});
sponIqr = cellfun(@(qs) diff(qs, 1, 1), sponQs, fnOpts{:});

respQs = cellfun(@(bs) quantile(bs(brFlag,:), [1,3]/4), behStack, fnOpts{:});
respMed = cellfun(@(bs) median(bs(brFlag,:), 1, "omitnan"), behStack, fnOpts{:});
respIqr = cellfun(@(bs) iqr(bs(brFlag,:), 1), behStack, fnOpts{:});

sponZ_bs = cellfun(@(x, xm, xs) myzs(x, xm, xs, 1), behStack, ...
    sponMu, sponSig, fnOpts{:});
allZ_bs = cellfun(@(bs) zscore(bs, 0, 1), behStack, fnOpts{:});

%%
pk_loc = cellfun(@(bs, m) getWaveformCriticalPoints(bs, fr), behStack, fnOpts{:});
pk_loc = cellfun(@(b) cellfun(@(t) t+bvWin(1), b, fnOpts{:}), pk_loc, fnOpts{:});
pk_rloc = cellfun(@(bs, m) getWaveformCriticalPoints(bs(brFlag,:), fr), behStack, fnOpts{:});
pk_rloc = cellfun(@(b) cellfun(@(t) t+brWin(1), b, fnOpts{:}), pk_rloc, fnOpts{:});
pk_sloc = cellfun(@(bs, m) getWaveformCriticalPoints(bs(bsFlag,:), fr), behStack, fnOpts{:});
pk_sloc = cellfun(@(b) cellfun(@(t) t+max(bvWin(1),bsWin(1)), b, fnOpts{:}), pk_sloc, fnOpts{:});
pvpt = cellfun(@(bs, pl, m) arrayfun(@(tr) interp1(behTx, bs(:,tr), ...
    pl{tr,1}, "cubic"), trigSubs, fnOpts{:}), behStack, pk_loc, fnOpts{:});
tend_mdl = cellfun(@(bs) arrayfun(@(tr) fit_poly(behTx, bs(:,tr), 1), ...
    trigSubs, fnOpts{:}), behStack, fnOpts{:});
right_side_flag = cellfun(@(pk) cellfun(@(tr) tr > 0, pk(:,1), ...
    fnOpts{:}), pk_loc, fnOpts{:});
frstPk = cellfun(@(b) cellfun(@(t) t(find(t>0.01, 1, "first")), b(:,1), ...
    fnOpts{:}), pk_loc, fnOpts{:});
frstPk_all = cellfun(@(b) cat(1, b{:}), frstPk, fnOpts{:});
[~, trOrd] = cellfun(@(b) sort(b), frstPk_all, fnOpts{:});
figure('Color', 'w'); hold on;
arrayfun(@(b) boxchart(b+zeros(size(frstPk_all{b})), frstPk_all{b}, ...
    bxOpts{:}, 'k'), 1:Nbs)
arrayfun(@(t) scatter(t + zeros(size(frstPk_all{t})) + ...
    random(jDist, size(frstPk_all{t})), frstPk_all{t}, 'k.'), 1:Nbs)
%%
bxFigs = gobjects(Nbs,1);
pdFlag = false(Ntr, Nbs); pmFlag = pdFlag; mvFlag = pdFlag;
pk_dist = cell(Ntr, Nbs); rosFlag = pk_dist; startleFlag = pdFlag;
froFlag = pdFlag; rgFlag = true(size(froFlag));
axs = gobjects(2,1); rTime = cell(Nbs, 1);
for cbs = 1:Nbs
    bxFigs(cbs) = figure("Name", behNames(cbs), fgOpts{:});
    axs(1) = subplot(10,1,1:8, "Parent", bxFigs(cbs));
    hold(axs(1), "on")
    for ctr = trigSubs(:)'
        if ~isempty(right_side_flag{cbs}{ctr}) && ...
                sum(right_side_flag{cbs}{ctr}) > 1
            sLoc = ctr+zeros(sum(~right_side_flag{cbs}{ctr}),1)-box_sep;
            sVal = pvpt{cbs}{ctr}(~right_side_flag{cbs}{ctr});% - ...
            %(pk_loc{cbs}{ctr,1}(~right_side_flag{cbs}{ctr}).^[1,0])*...
            %tend_mdl{cbs}{ctr};
            rLoc = ctr+zeros(sum(right_side_flag{cbs}{ctr}),1)+box_sep;
            rVal = pvpt{cbs}{ctr}(right_side_flag{cbs}{ctr});% - ...
            %(pk_loc{cbs}{ctr,1}(right_side_flag{cbs}{ctr}).^[1,0])*...
            %tend_mdl{cbs}{ctr};
            rTime = pk_loc{cbs}{ctr,1}(right_side_flag{cbs}{ctr});
            if ~isempty(sVal) && ~isempty(rVal)
                [~, pmFlag(ctr, cbs)] = ranksum(sVal, rVal);
                sIqr = iqr(sVal); sQs = quantile(sVal,[1,3]/4);
                rosFlag{ctr, cbs} = isWhiskOutlier(rVal, sQs(:), sIqr);
                frPkSubs = find(rTime>0.01,3,"first");
                froFlag(ctr, cbs) = any(rosFlag{ctr,cbs}(frPkSubs));
                pk_dist{ctr, cbs} = distmatrix(sVal, rVal);
                if numel(sVal) > 1 && numel(rVal) > 1
                    pdFlag(ctr, cbs) = ansaribradley(sVal-median(sVal), ...
                        rVal-median(rVal));
                    if froFlag(ctr, cbs)
                        % Is outlier at least a 10th of the signal range?
                        rgFlag(ctr, cbs) = mean(abs((cumsum(any((...
                            pk_dist{ctr, cbs}(:, frPkSubs) > (mvRng(cbs)/10)), ...
                            2)).^[1,0]) * fit_poly([0, size(sVal,1)], [0,1], 1))) ...
                            > 0.5;
                    end
                end
            end
            boxchart(axs(1), sLoc, sVal, bxOpts{:}, 'g')
            boxchart(axs(1), rLoc, rVal, bxOpts{:}, 'r')
            scatter(axs(1), sLoc+random(jDist,size(sLoc)), sVal,'g.')
            scatter(axs(1), rLoc+random(jDist,size(rLoc)), rVal,'r.')
        end
    end
    mvFlag(:, cbs) = (pdFlag(:,cbs) | pmFlag(:,cbs) | froFlag(:,cbs)) & ...
        rgFlag(:, cbs);
    axs(2) = subplot(10,1,9:10); stem(axs(2), mvFlag(:,cbs))
    xticks(axs, trigSubs); xticklabels(axs, trigSubs); set(axs, axOpts{:})
    linkaxes(axs, 'x')
    figure('Color', 'w'); hold on; arrayfun(@(tr) patch([behTx; NaN], ...
        [behStack{cbs}(:,tr) - sponMed{cbs}(tr); NaN], ...
        [repmat(tr,Nbt,1);NaN], [ones(Nbt,1);nan], ...
        'EdgeColor', 'r', 'EdgeAlpha', 0.5), find(mvFlag(:,cbs)))
    arrayfun(@(tr) patch([behTx; NaN], ...
        [behStack{cbs}(:,tr) - sponMed{cbs}(tr); NaN], ...
        [repmat(tr,Nbt,1);NaN], [ones(Nbt,1);nan], ...
        'EdgeColor', 'k', 'EdgeAlpha', 0.3), find(~mvFlag(:,cbs)))
    text(repmat(behTx(1),1,Ntr), zeros(Ntr,1), 1:Ntr, string((1:Ntr)'), ...
        'HorizontalAlignment','right')
    patch([0,0,0,0], [min(ylim), min(ylim), max(ylim), max(ylim)], ...
        [0, Ntr, Ntr, 0], [1,1,1,1], 'EdgeColor', 'none', ...
        'FaceAlpha', 0.2, 'FaceColor', 'k')
    %{
    temp_x = behStack{cbs}(bsFlag | brFlag,:);
    figure; ; boxplot(temp_x(:), bxGrup(:), bxOpts{:});
    title(behNames(cbs)); hold on; scatter(2:2:2*Ntr,mvpt{cbs},'rx');
    xticks(2:2:2*Ntr); xticklabels(1:Ntr)
    xsFlag(:, cbs) = isWhiskOutlier(mvpt{cbs}, sponQs{cbs}, sponIqr{cbs});
    xrFlag(:, cbs) = isWhiskOutlier(mvpt{cbs}, respQs{cbs}, respIqr{cbs});
    ssFlag(:, cbs) = ansaribradley(behStack{cbs}(bsFlag,:)-sponMed{cbs}, ...
        behStack{cbs}(brFlag,:)-median(behStack{cbs}(brFlag,:),1));
    qr_mi(:, cbs) = getMI(sponIqr{cbs}, respIqr{cbs});
    mxs_mi(:, cbs) = getMI(sponMed{cbs}(:), mvpt{cbs});
    mxr_mi(:, cbs) = getMI(respMed{cbs}(:), mvpt{cbs});
    %zst(:, cbs) = xwFlag(:,cbs) & ssFlag(:,cbs);
    zst(:, cbs) = ...
        (((msFlag(:,cbs) | ssFlag(:,cbs)) & abs(qr_mi(:, cbs)) > 0.1) | ... State change
        ((xsFlag(:,cbs) & xrFlag(:, cbs)) & abs(mxs_mi(:, cbs)) > 0.1)) & ... Short peak
        cellfun(@(x) ~isempty(x), pk_rloc{cbs}(:,1));    %
    subplot(10,1,9:10); stem(2:2:2*Ntr, abs(zst(:, cbs)))
    xticks(2:2:2*Ntr); xticklabels(1:Ntr); linkaxes(get(gcf, 'Children'), 'x')
    figure; hold on; arrayfun(@(tr) patch([behTx; NaN], ...
        [behStack{cbs}(:,tr) - sponMed{cbs}(tr); NaN], ...
        [repmat(tr,Nbt,1);NaN], [ones(Nbt,1);nan], ...
        'EdgeColor', 'r', 'EdgeAlpha', 0.5), find(zst(:,cbs)))
    arrayfun(@(tr) patch([behTx; NaN], ...
        [behStack{cbs}(:,tr) - sponMed{cbs}(tr); NaN], ...
        [repmat(tr,Nbt,1);NaN], [ones(Nbt,1);nan], ...
        'EdgeColor', 'k', 'EdgeAlpha', 0.3), find(~zst(:,cbs)))
    text(repmat(behTx(1),1,Ntr), zeros(Ntr,1), 1:Ntr, string((1:Ntr)'), ...
        'HorizontalAlignment','right')
    % xline(0,'k:')
    patch([0,0,0,0], [min(ylim),min(ylim),max(ylim),max(ylim)], ...
        [0, Ntr, Ntr, 0], [1,1,1,1], 'EdgeColor', 'none', ...
        'FaceAlpha', 0.2, 'FaceColor', 'k')
    end
    %}
end

%%
for cbs = 1:Nbs
    for ctr = 1:Ntr
        figure; plot(behTx, behStack{cbs}(:,ctr), 'k');
        title(behNames(cbs)+" trial "+string(ctr))
        hold on;
        scatter(pk_rloc{cbs}{ctr,1}, interp1(behTx, behStack{cbs}(:,ctr), ...
            pk_rloc{cbs}{ctr,1}, "cubic"), 'ro')
        scatter(pk_rloc{cbs}{ctr,2}, interp1(behTx, behStack{cbs}(:,ctr), ...
            pk_rloc{cbs}{ctr,2}, "cubic"), 'r.')
        scatter(pk_sloc{cbs}{ctr,1}, interp1(behTx, behStack{cbs}(:,ctr), ...
            pk_sloc{cbs}{ctr,1}, "cubic"), 'go')
        scatter(pk_sloc{cbs}{ctr,2}, interp1(behTx, behStack{cbs}(:,ctr), ...
            pk_sloc{cbs}{ctr,2}, "cubic"), 'g.')
        scatter(pk_loc{cbs}{ctr,1}, interp1(behTx, behStack{cbs}(:,ctr), ...
            pk_loc{cbs}{ctr,1}, "cubic"), 10, 'ko')
        scatter(pk_loc{cbs}{ctr,2}, interp1(behTx, behStack{cbs}(:,ctr), ...
            pk_loc{cbs}{ctr,2}, "cubic"), 'k.')
        xline(0, 'k:')
    end
end



%% Z-score measurements
% Z-score w.r.t.
jDist = makedist('Normal', 'mu', 0, 'sigma', 1/8);

%{
[mvpt, max_time] = cellfun(@(bs) getMaxAbsPerTrial(bs, brWin, behTx), behStack, ...
    fnOpts{:});
mvpt_tz = cellfun(@(bs) getMaxAbsPerTrial(bs, brWin, behTx), sponZ_bs, ...
    fnOpts{:});
mvpt_az = cellfun(@(bs) getMaxAbsPerTrial(bs, brWin, behTx), allZ_bs, ...
    fnOpts{:});
%}
figure; boxchart([max_time{:}], 'MarkerStyle', 'none');
hold on; scatter(repmat(1:Nbs, Ntr, 1) + random(jDist,[Ntr, Nbs]), ...
    [max_time{:}], 'kx')
% Whisker outlier: is the point ssmaller than Q1 - 1.5*IQR or greater than
% Q3 + 1.5*IQR?


bxGrup = reshape(1:2*Ntr,2,[]);
bxGrup = [repmat(bxGrup(1,:),sum(bsFlag), 1); repmat(bxGrup(2,:),sum(brFlag), 1)];
% Normality assumption
szFlag = cellfun(@(z) cell2mat(arrayfun(@(a) abs(z)>a, signTh, fnOpts{:})), ...
    mvpt_tz, fnOpts{:});
azFlag = cellfun(@(z) cell2mat(arrayfun(@(a) abs(z)>a, signTh, fnOpts{:})), ...
    mvpt_az, fnOpts{:});
% Non-parametric tests
zst = false(Ntr, Nbs); xsFlag = false(Ntr, Nbs); ssFlag = xsFlag;
xrFlag = xsFlag;
[~, msFlag] = cellfun(@(bs) arrayfun(@(tr) ...
    ranksum(bs(bsFlag,tr), bs(brFlag,tr)), (1:Ntr)'), behStack, fnOpts{:});
msFlag = [msFlag{:}]; qr_mi = zeros(Ntr, Nbs); mxs_mi = qr_mi; mxr_mi = qr_mi;
bxOpts = {'PlotStyle', 'compact', 'Notch', 'on'};

for cbs = 1:Nbs
    temp_x = behStack{cbs}(or(bsFlag, brFlag),:);
    figure; subplot(10,1,1:8); boxplot(temp_x(:), bxGrup(:), bxOpts{:});
    title(behNames(cbs)); hold on; scatter(2:2:2*Ntr,mvpt{cbs},'rx');
    xticks(2:2:2*Ntr); xticklabels(1:Ntr)
    xsFlag(:, cbs) = isWhiskOutlier(mvpt{cbs}, sponQs{cbs}, sponIqr{cbs});
    xrFlag(:, cbs) = isWhiskOutlier(mvpt{cbs}, respQs{cbs}, respIqr{cbs});
    ssFlag(:, cbs) = ansaribradley(behStack{cbs}(bsFlag,:)-sponMed{cbs}, ...
        behStack{cbs}(brFlag,:)-median(behStack{cbs}(brFlag,:),1));
    qr_mi(:, cbs) = getMI(sponIqr{cbs}, respIqr{cbs});
    mxs_mi(:, cbs) = getMI(sponMed{cbs}(:), mvpt{cbs});
    mxr_mi(:, cbs) = getMI(respMed{cbs}(:), mvpt{cbs});
    %zst(:, cbs) = xwFlag(:,cbs) & ssFlag(:,cbs);
    zst(:, cbs) = ...
        (((msFlag(:,cbs) | ssFlag(:,cbs)) & abs(qr_mi(:, cbs)) > 0.1) | ... State change
        ((xsFlag(:,cbs) & xrFlag(:, cbs)) & abs(mxs_mi(:, cbs)) > 0.1)) & ... Short peak
        cellfun(@(x) ~isempty(x), pk_rloc{cbs}(:,1));    %
    subplot(10,1,9:10); stem(2:2:2*Ntr, abs(zst(:, cbs)))
    xticks(2:2:2*Ntr); xticklabels(1:Ntr); linkaxes(get(gcf, 'Children'), 'x')
    figure; hold on; arrayfun(@(tr) patch([behTx; NaN], ...
        [behStack{cbs}(:,tr) - sponMed{cbs}(tr); NaN], ...
        [repmat(tr,Nbt,1);NaN], [ones(Nbt,1);nan], ...
        'EdgeColor', 'r', 'EdgeAlpha', 0.5), find(zst(:,cbs)))
    arrayfun(@(tr) patch([behTx; NaN], ...
        [behStack{cbs}(:,tr) - sponMed{cbs}(tr); NaN], ...
        [repmat(tr,Nbt,1);NaN], [ones(Nbt,1);nan], ...
        'EdgeColor', 'k', 'EdgeAlpha', 0.3), find(~zst(:,cbs)))
    text(repmat(behTx(1),1,Ntr), zeros(Ntr,1), 1:Ntr, string((1:Ntr)'), ...
        'HorizontalAlignment','right')
    % xline(0,'k:')
    patch([0,0,0,0], [min(ylim),min(ylim),max(ylim),max(ylim)], ...
        [0, Ntr, Ntr, 0], [1,1,1,1], 'EdgeColor', 'none', ...
        'FaceAlpha', 0.2, 'FaceColor', 'k')
end
%% Indexing p(r), p(s), and p(r,s)
itR = false(size(behTx_all)); itS = false(size(itR));
atR = false(size(behTx_all)); atS = false(size(atR));
trial_duration_samples = round(diff(brWin)*fr);
trial_offset_samples = round(brWin(1)*fr);

itSubs = round(itTimes{iSub}(:,1)*fr);
atSubs = round(atTimes{aSub}*fr);
% Intan subscripts
% Response window
itRespDur_subs = arrayfun(@(tr) itSubs(tr) + ...
    (0:trial_duration_samples)' + trial_offset_samples, trigSubs, fnOpts{:});
itRespDur_subs = cat(1, itRespDur_subs{:});
itR(itRespDur_subs) = true;
% Spontaneous window
itSponDur_subs = arrayfun(@(tr) itSubs(tr) + ...
    (-trial_duration_samples:0)' - trial_offset_samples, trigSubs, fnOpts{:});
itSponDur_subs = cat(1, itSponDur_subs{:});
itS(itSponDur_subs) = true;
% Arduino subscripts
% Response window
atRespDur_subs = arrayfun(@(tr) atSubs(tr) + ...
    (0:trial_duration_samples)' + trial_offset_samples, trigSubs, fnOpts{:});
atRespDur_subs = cat(1, atRespDur_subs{:});
atR(atRespDur_subs) = true;
% Spontaneous window
atSponDur_subs = arrayfun(@(tr) atSubs(tr) + ...
    (-trial_duration_samples:0)' - trial_offset_samples, trigSubs, fnOpts{:});
atSponDur_subs = cat(1, atSponDur_subs{:});
atS(atSponDur_subs) = true;


%% Probabilities
rngDLC = range(behDLCSignals, 1); Nbns = 256;
mn = min(behDLCSignals, [], 1); mx = max(behDLCSignals, [], 1);
rngRS = range(vf, "all"); mnr = min(vf); mxr = max(vf);
lmt = mat2cell([mn(:),mx(:)], ones(size(mn,2),1));
hsOpts = {'Normalization', 'pdf', 'BinWidth', 'BinLimits'};
hs2Opts = [repmat({'XBinLimits'},size(mn,2),1), lmt, ...
    repmat({'YBinLimits'},size(mn,2),1), lmt];

[p_r, d_r] = arrayfun(@(s) histcounts(behDLCSignals(itRespDur_subs,s), ...
    hsOpts{1:3}, rngDLC(s)/Nbns, hsOpts{4}, [mn(s), mx(s)]), ...
    1:size(behDLCSignals,2), fnOpts{:});
[p_rr, d_rr] = histcounts(vf(atRespDur_subs),hsOpts{1:3}, rngRS/Nbns, ...
    hsOpts{4}, [mnr, mxr]);
[p_s, d_s] = arrayfun(@(s) histcounts(behDLCSignals(itSponDur_subs,s), ...
    hsOpts{1:3}, rngDLC(s)/Nbns, hsOpts{4}, [mn(s), mx(s)]), ...
    1:size(behDLCSignals,2), fnOpts{:});
[p_sr, d_sr] = histcounts(vf(atSponDur_subs), hsOpts{1:3}, rngRS/Nbns, ...
    hsOpts{4}, [mnr, mxr]);
[p_rs, d_rs] = arrayfun(@(s) histcounts2(behDLCSignals(itRespDur_subs,s), ...
    behDLCSignals(itSponDur_subs,s), hsOpts{1:3}, repmat(rngDLC(s)/Nbns,2,1), ...
    hs2Opts{s,:}), 1:size(behDLCSignals,2), fnOpts{:});
p_r = [p_r(:)', {p_rr}]; p_s = [p_s(:)', {p_sr}];
kld = cellfun(@(a,b) KullbackLeiblerDivergence(a, b), p_r, p_s);
minfo = arrayfun(@(s) mi_cont_cont(behDLCSignals(itSponDur_subs,s), ...
    behDLCSignals(itRespDur_subs,s),1), 1:size(mn,2));

