fr = VideoReader(fullfile(behDir, "roller2021-06-23T16_16_22.avi"));
fr = fr.FrameRate;
rp = readRollerPositionsFile(fullfile(behDir, ...
    "Roller_position2021-06-23T16_16_07.csv"));
[vf, rollTx] = getRollerSpeed(rp, fr);
pTms = getCSVTriggers(fullfile(behDir, "Puff.csv")) - rollTx(1);
lTms = getCSVTriggers(fullfile(behDir, "Laser.csv")) - rollTx(1);
trialID = readCSVtimeStamps(fullfile(behDir, "PL.csv"));
trialFlag = false(size(pTms));
trialFlag(sub2ind(size(pTms), 1:size(pTms,1), ...
    trialID(any(trialID(:,1) == [1,2], 2), 1)')) = true;

timeLapse = [-0.25, 0.5];

[~, vStack] = getStacks(false, round(pTms*fr), 'on', timeLapse, fr, fr, ...
    [], vf);
[~, Nt, Na] = size(vStack);
stMdl = fit_poly([1, Nt], timeLapse, 1); stTx = ((1:Nt)'.^[1,0]) * stMdl;
figure; plot(stTx, squeeze(vStack(:,:,trialFlag(:,1))))