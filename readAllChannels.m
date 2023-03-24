afOpt = {'UniformOutput', false};

%% Reading the binary file
clWaveforms = cell(1,2);
% Taking ~1.25 ms around the spike.
spikeWaveTime = 2*round(1.25e-3 * fs) + 1;
spikeSamples = (spikeWaveTime - 1)/2;

% pcFeat = readNPY(fullfile(dataDir, 'pc_features.npy'));
% pcInd = readNPY(fullfile(dataDir, 'pc_feature_ind.npy'));
spkSubs = round(spkTms.*fs);
% [chs2read, readOrder, repeatChs] = unique(ch2read);
chs2read = 20; ch2read = 20;
answ = 1;

cchan = 1;
% Main loop
while ~feof(fID) && cchan <= size(chs2read,1)
    % Computing the location of the channel features
    % pcIdx = ch2read(cchan) == chanMap(pcInd(clTempSubs{cchan}+1,:)+1);
    % clFeat = pcFeat(spkIdx(:,cchan), :, pcIdx);
    
    % Jumping to the channel
    
    % Getting all clusters from the considered channel
    clustChanIdx = ch2read == chs2read(cchan); Nccl = sum(clustChanIdx);
    Nspks = size(spkSubs, 1);
    spkLbls = arrayfun(@(x) x*ones(Nspks(x),1), 1:Nccl,...
        afOpt{:}); spkLbls = cat(1, spkLbls{:});
    chSpks = [cat(1, spkSubs), spkLbls];
    [ordSpks, spkOrd] = sort(chSpks(:,1), 'ascend');
    % Computing the distance from spike to spike
    spkDists = [ordSpks(1);diff(ordSpks)];
    % Allocating space for the spikes
    waveform = zeros(spikeWaveTime, 64, sum(Nspks), 'single');
    %fig = figure('Color',[1,1,1],'Visible', 'off');
    %ax = axes('Parent', fig); ax.NextPlot = 'add';
    %subSet = 1:floor(numel(spkDists)*0.1);
    for cspk = 1:sum(Nspks)
        % Jumping to 1 ms before the time when the spike occured
        fseek(fID, 2*((Nch+1)*(spkDists(cspk) - spikeSamples)), 'cof');
        % Reading the waveform
        cwf = fread(fID, [spikeWaveTime, 64], 'int16=>single');
        % Assigning the waveform to the saving variable. If the waveform is
        % cut, then assign only the gathered piece.
        waveform(:,:,cspk) = cwf;
        % Jumping back to the exact time of the spike
        fseek(fID, -2*((Nch+1)*(spikeSamples+1)), 'cof');
        %    if ismember(cspk,subSet)
        %        plot(ax,waveform(:,cspk),'DisplayName',num2str(cspk));
        %    end
    end
    fprintf(1,' done!\n')
    clWaveforms(clustChanIdx,:) = [clusterID(clustChanIdx),...
        arrayfun(@(x) waveform(:,spkOrd(chSpks(:,2) == x)), (1:Nccl)', afOpt{:})];
    cchan = cchan + 1;
    frewind(fID);
    %fig.Visible = 'on';
end
fclose(fID);
