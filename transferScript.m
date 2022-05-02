hold on; plot((0:size(vx,2)-1)/fr + rollTx(1), vx)
figure; plot((0:size(vx,2)-1)/fr + rollTx(1), [vx', vx_f'])
fclose(fID);
[~] = fclose(fID);
figure; cwt(vf, 'amor', fr)
[b9, a9] = butter(5, (2*18)/fr, 'low');
vf = filtfilt(b,a,vx);
vx_f2 = filtfilt(b,a,double(vx));
figure; plot((0:size(vx,2)-1)/fr + rollTx(1), [vx', vx_f', vx_f2'])
figure; cwt(vx_f2, 'amor', fr)
[~] = fclose(fID);
fID = fopen("Z:\Emilio\SuperiorColliculusExperiments\Roller\Batch1_ephys\GAD48_S1\210625\Behaviour\Laser.csv", 'r');
lns = textscan(fID,'%s',Inf);
[~] = fclose(fID);
lTms = cellfun(@(x) datetime(cleanStr(x, ','), 'InputFormat', 'uuuu-MM-dd''T''HH:mm:ss.SSS'), lns{:});
fID = fopen("Z:\Emilio\SuperiorColliculusExperiments\Roller\Batch1_ephys\GAD48_S1\210625\Behaviour\PL.csv",'r');
lns = textscan(fID,'%s',Inf); lns = lns{:}; trialID = lns;
[~] = fclose(fID);
trialSub = str2double(extractBefore(trialID, ','));