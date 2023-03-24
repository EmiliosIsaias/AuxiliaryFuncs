%%
nDist = makedist('Normal', "mu", 0, "sigma", 1);
alph = 1 - [0.2, 0.1 ,0.05,0.001];
signTh = arrayfun(@(z) fminbnd(@(y) ...
    norm(integral(@(x) pdf(nDist, x), -y, y) - z, 2), -6, 6), alph);
myzs = @(x, mu, st) (x - mu)./(st.*(st~=0) + 1.*(st==0));
sponMu = cellfun(@(bs) mean(bs(bsFlag,:), 1, "omitnan"), behStack, fnOpts{:});
sponSig = cellfun(@(bs) std(bs(bsFlag,:), 0, 1, "omitnan"), behStack, fnOpts{:});
respSig = cellfun(@(bs) std(bs(brFlag,:), 0, 1, "omitnan"), behStack, fnOpts{:});
trialZ_bs = cellfun(myzs, behStack, sponMu, sponSig, fnOpts{:});

allZ_bs = cellfun(@(bs) zscore(bs, 0, 1), behStack, fnOpts{:});


%% Z-score measurements
% Z-score w.r.t.
mvpt = cellfun(@(bs) getMaxAbsPerTrial(bs, brWin, behTx), behStack, ...
    fnOpts{:});
mvpt_tz = cellfun(@(bs) getMaxAbsPerTrial(bs, brWin, behTx), trialZ_bs, ...
    fnOpts{:});
mvpt_az = cellfun(@(bs) getMaxAbsPerTrial(bs, brWin, behTx), allZ_bs, ...
    fnOpts{:});

zFlag = cellfun(@(z) abs(z)>1.96, mvpt_tz, fnOpts{:});
trialSig = cellfun(@(bs) std(bs-median(bs(bsFlag,:),1), 0, 2), behSgnls, ...
    fnOpts{:});

szFlag = cellfun(@(z) cell2mat(arrayfun(@(a) abs(z)>a, signTh, fnOpts{:})), ...
    mvpt_tz, fnOpts{:});
azFlag = cellfun(@(z) cell2mat(arrayfun(@(a) abs(z)>a, signTh, fnOpts{:})), ...
    mvpt_az, fnOpts{:});
for cbs = 1:Nbs
    figure; boxplot([behStack{cbs}(bsFlag,:); mvpt{cbs}(:)']);
    title(behNames(cbs)); hold on; scatter(1:Ntr,mvpt{cbs},'go')
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

