pairedStim = arrayfun(@(x) Conditions(1).Triggers(:,1) == ...
Conditions(x).Triggers(:,1)', [3,6,7], fnOpts{:});
pairedStim = cellfun(@(x) any(x, 2), pairedStim, fnOpts{:});
pairedStim = cat(2, pairedStim{:});

timeLapse = [-50, 300] * 1e-3;
chNames = [13,26,40,51];
shNames = 1:4;
stTx = (1:size(LFPstack,2))'/fs - 0.05;
for cch = 1:size(data_scaled,2)
    for ccond = 1:2
        figure; plot(stTx, squeeze(LFPstack(cch,:,pairedStim(:,ccond))), "LineWidth", 0.005, "Color", 0.75*ones(1,3))
        hold on; plot(stTx, mean(squeeze(LFPstack(cch,:,pairedStim(:,ccond))),2), "LineWidth", 1, "Color", "k")
        yyaxis right; plot(stTx, mean(squeeze(lsStack(1,:,pairedStim(:,ccond))),2), "Color", 'b')
        title(sprintf('Channel %d, shank %d', chNames(cch), shNames(cch)))
        set(gca, "Box", "off", "Color", "none")
        xlim(timeLapse)
    end
end

trMvFlag = arrayfun(@(cr) behRes(1).Results(cr).MovStrucure.MovmentFlags, ...
    1:size(behRes(1).Results,2), fnOpts{:}); trMvFlag = cat(3, trMvFlag{:});
BIscale = sum(trMvFlag,3);
BIscale = arrayfun(@(cc) BIscale(pairedStim(:,cc), cc), 1:Nccond, ...
    fnOpts{:});
[hg, hg_bin] = cellfun(@(c) histcounts(c, "BinMethod", "integers"), ...
    BIscale, fnOpts{:}); hg = cat(1, hg{:}); hg_bin = cat(1, hg_bin{:});

%%
fnOpts = {'UniformOutput', false};
getChildFolder = @(x) fullfile(x.folder, x.name);
miceDir = dir(fullfile("Z:\Emilio\SuperiorColliculusExperiments\Roller\Batch14_ephys.MC", "*\WT*"));
expDirs = arrayfun(@(d) dir(fullfile(getChildFolder(d), '23*')), miceDir, fnOpts{:});
for cm = 1:numel(expDirs)
    mBehRes = analyseBehaviour_allSessions(expDirs{cm});
    consCondNames = cellfun(@(ms) {ms(:).ConditionName}, mBehRes, fnOpts{:});
    movProp_perSess = cellfun(@(ms) arrayfun(@(bp) ...
        ms(1).Results(bp).MovStrucure.MovmentFlags, 1:4, fnOpts{:}), ...
        mBehRes, fnOpts{:});
    movProp_perSess = cellfun(@(mps) cat(3, mps{:}), movProp_perSess, fnOpts{:});
    mvCnt = cellfun(@(x) squeeze(sum(x, 1)), movProp_perSess, fnOpts{:});
    
    Ntc = cellfun(@(mc) [mc(:).NTrials], mBehRes, fnOpts{:});

    mvCnt2 = cellfun(@(cnt) cnt([1,2],:), mvCnt, fnOpts{:});
    mvCnt2 = cat(3, mvCnt2{:});

    Ntc2 = cellfun(@(tot) tot([1,2]), Ntc, fnOpts{:});
    Ntc2 = cat(1, Ntc2{:});

    mouseMovProp = sum(mvCnt2,3)./sum(Ntc2)';

   
    behSigName = {'Stimulated whiskers','Non-stimulated whiskers', ...
        'Nose', 'Roller speed'};
    tmpBehRes = [];
    for cc = 1:2
        tmpBehRes = [tmpBehRes, struct('ConditionName', consCondNames{1}{cc}, ...
            'Results', struct('BehSigName',[],'MovProbability',0))];
        for cbs = 1:4
            auxStruct = struct('BehSigName', behSigName{cbs}, ...
                'MovProbability', mouseMovProp(cc, cbs));
            tmpBehRes(cc).Results(cbs) = auxStruct;
        end
    end
    [pAreas, ~, polyFig] = createBehaviourIndex(tmpBehRes);
    biFigPttrn = sprintf("%s BehIndex%%s", miceDir(cm).name);
    biFigPttrn = sprintf(biFigPttrn, sprintf(" %s (%%.3f)", consCondNames{cm}{[1,2]}));
    tmpBehRes = arrayfun(@(bs, ba) setfield(bs,'BehIndex', ba), tmpBehRes, pAreas);
    set(polyFig, 'UserData', {tmpBehRes, mBehRes})
    biFN = sprintf(biFigPttrn, pAreas);
    saveFigure(polyFig, fullfile(expDirs{cm}(1).folder, biFN), true, true);
    save(fullfile(expDirs{cm}(1).folder, 'AllSessions.mat'), 'mBehRes', 'tmpBehRes')
end