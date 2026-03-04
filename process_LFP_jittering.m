%% Parameters
fnOpts = {'UniformOutput', false};
ccvv = @(x) cat(1,x{:});
vw = [-1,5]*1e-2;
bin_size = 1e-4;
n_bins = round(diff(vw)/bin_size);
tx_psth = ((1:n_bins)-0.5)*bin_size + vw(1);

fullpath = @(x) fullfile(x.folder, x.name);
fs_low = 500;
[b,a] = butter(2,2*[0.2,20]/fs_low,'bandpass');
jittering_folder = fullfile('Z:','Jesus','Jittering','FULLYCurated');
jittering_selected_experiments = [3:2:13,14:2:18,22,27];
all_jittering_exps = dir(jittering_folder);
all_jittering_exps(~[all_jittering_exps.isdir]) = [];
non_exp_folders = cellfun(@(x) isempty(x), ...
    regexp({all_jittering_exps.name},'^\d{2}_'));
all_jittering_exps(non_exp_folders) = [];
experiments_lfp = all_jittering_exps(jittering_selected_experiments);
%% Loop
for ce = experiments_lfp'
    % Load necessary files
    trig_file = dir(fullfile(fullpath(ce),'**','*analysis.mat'));
    if numel(trig_file) > 1
        [~,bigger_file] = max([trig_file.bytes]);
        trig_file = trig_file(bigger_file);
    end
    load(fullpath(trig_file), "Conditions", "Triggers");
    fs_file = dir(fullfile(fullpath(ce),'**', '*_sampling_frequency.mat'));
    if numel(fs_file)>1
        fs_file = fs_file(1);
    end
    load(fullpath(fs_file),'fs')
    spk_file = dir(fullfile(fullpath(ce),"**","*_Spike_Times.mat"));
    if numel(spk_file) > 1
        spk_file = spk_file(1);
    end
    load(fullpath(spk_file))
    % Downsample and filter the LFP
    lfp = Triggers.LFP; clearvars Triggers;
    lfp_ds = resample(lfp,(((1:length(lfp))-0.5)/fs),fs_low);
    tx = ((1:length(lfp_ds))-0.5)/fs_low;
    lfp = filtfilt(b,a,lfp_ds);
    temp = zscore(lfp);
    brain_state = fitSpline(tx, temp, 1, 0.65, 0.75, 1);
    nf_flag = brain_state <= 0;
    % Assaign a state per trial
    sObj = StepWaveform(nf_flag,fs_low);
    sObj.MinIEI = 0.5;
    nfs_triggers = sObj.Triggers/fs_low;
    wscc = contains({Conditions.name}, {'whisker','piezo'}, ...
        "IgnoreCase", true) & contains({Conditions.name}, 'Control');
    ws_triggers = sort(cat(1, Conditions(wscc).Triggers), "ascend")/fs;
    ws_triggers = ws_triggers(:,1);
    on_stim_flag = any(ws_triggers' > nfs_triggers(:,1) &...
        ws_triggers' < nfs_triggers(:,2), 1);
    if all(~on_stim_flag)
        continue;
    end
    % Compute for how long has the LFP been in 'upstate' (1) or 'downstate'
    % (0)
    t1_ids = find(on_stim_flag);
    t2_ids = find(~on_stim_flag);
    delta_t1 = ws_triggers(:,1)' - nfs_triggers(:,1);
    [~, on_sub] = min(abs(delta_t1));
    temp = arrayfun(@(x) delta_t1(on_sub(x),x), 1:numel(on_sub));
    [t1_time, t1_ord] = sort(temp(on_stim_flag), "descend");
    % Special treatment for off boundaries. Only considering off 'upstates'
    % that are before trials.
    delta_t2 = ws_triggers(:,1)' - nfs_triggers(:,2);
    delta_t2(delta_t2<0) = nan;
    [~, off_sub] = min(delta_t2);
    temp = arrayfun(@(x) delta_t2(off_sub(x),x), 1:numel(off_sub));
    [t2_time, t2_ord] = sort(temp(~on_stim_flag), "ascend");
    % Compute PSTHs and first spike per unit, per trial
    PSTH = zeros(numel(spike_times), ... Number of units
        n_bins,...Number of bins for the given viewing window and bin size
        size(ws_triggers,1), ...         Number of trials
        "single");
    fst = zeros(numel(spike_times),size(ws_triggers,1));
    for ct = 1:size(ws_triggers,1)
        % Get PSTH per trial and unit
        temp = cellfun(@(x) histcounts(x-ws_triggers(ct,1), ...
            'BinLimits', vw, 'BinWidth',bin_size), ...
            spike_times, fnOpts{:});
        temp = ccvv(temp);
        PSTH(:,:,ct) = temp;
        % Get relative spike times w.r.t. each trial per unit
        temp = cellfun(@(x) x - ws_triggers(ct,1), spike_times, fnOpts{:});
        temp = cellfun(@(x) min(x(x>2e-3 & x<5e-2)), temp, fnOpts{:});
        temp2 = cellfun(@(x) ~isempty(x), temp);
        fst(temp2,ct) = [temp{:}];
    end
    fst_ut = fst(:,[t1_ids(t1_ord),t2_ids(t2_ord)]);
    fst_ut(fst_ut==0) = nan;
    fst_median = median(fst_ut,1,"omitmissing");
    ns_flag = isnan(fst_median);
    [cw_spect, fax] = cwt(lfp_ds, fs_low);

    %% Plotting results
    f = figure("Color","w"); t = createtiles(f,1,1);
    ax = nexttile(t);
    imagesc(ax,tx_psth*1e3,[],squeeze(mean(PSTH(:,:,[t1_ids(t1_ord),t2_ids(t2_ord)]),1))')
    yticks(ax,1:numel(on_stim_flag))
    yticklabels(ax,[t1_time,t2_time])
    colormap(ax,inferno); cleanAxis(ax);
    set(ax,"TickDir","out"); xlabel(ax,"Time [ms]"); 
    ylabel(ax,'Down \leftarrow{} Up-state ''delay'' [s] \rightarrow{} Up', ...
        'Interpreter','tex')
    yline(ax,max(t1_ord)+0.5, 'w-')
    title(t,'Trial-by-trial population response')
    saveFigure(f,fullfile(spk_file.folder,'Figures', ...
        'Trial-by-trial responses'),1,0)
    f = figure("Color","w"); t = createtiles(f,1,1);
    ax = nexttile(t);
    histogram2(ax,repmat([t1_time,-t2_time]',numel(spike_times),1),fst_ut(:), ...
        [64,64], "FaceColor","flat","EdgeColor","none", "ShowEmptyBins","on"); 
    colormap(ax,inferno); hold(ax,"on")
    scatter(ax,[t1_time,-t2_time], fst_ut, 'w.')
    view(ax,0,90)
    set(ax,'TickDir','out'); xlabel('Up-sate "delay" [s]');
    xlabel(ax,'Up-state distribution [s]')
    ylabel(ax,'First spike time [s]')
    title(t,'First spike vs up-state "delay" distribution')
    cleanAxis(ax);
    saveFigure(f,fullfile(spk_file.folder,'Figures', ...
        'FST vs up-state distribution'),1,0)
    f = figure("Color","w"); t = createtiles(f,4,1);
    axs(1) = nexttile(t,1,[3,1]);
    imagesc(axs(1),tx,log10(fax),20*log10(abs(cw_spect)));
    yticklabels(axs(1), 10.^yticks(axs(1)));
    colormap(axs(1),inferno); ylabel(axs(1),'Frequency [Hz]')
    yyaxis(axs(1),'right'); line(axs(1),tx,lfp_ds,'Color','k'); 
    ylim(axs(1), mean(lfp_ds) + [-5,5]*std(lfp_ds))
    set(axs(1).XAxis, "Visible", "off")
    axs(2) = nexttile(t,4);
    histogram(axs(2),cat(1,spike_times{:}),((0:round(Nt/0.01)))*0.01, ...
        "FaceColor","k","EdgeColor","none")
    yyaxis(axs(2),"right"); 
    stem(axs(2), ws_triggers([t1_ids(t1_ord),t2_ids(t2_ord)],1), fst_median)
    y = axs(2).YLim;
    y2 = repmat(y(:), [1,size(nfs_triggers,1)]);
    y2 = padarray(y2,2,"symmetric","post");
    x = padarray(nfs_triggers,[0,1],"symmetric","both")';
    patch(axs(2),x,y2,ones(size(x)), 'EdgeColor', 'none', ...
        'FaceColor', 'b', 'FaceAlpha', 0.4)
    set(axs,"TickDir","out")
    yticklabels(axs(1),10.^yticks(axs(1)))
    yyaxis(axs(2),'left'); ylabel(axs(2),'MUA spike count')
    yyaxis(axs(2),'right'); ylabel(axs(2),'Median first spike')
    cb = colorbar(axs(1),"Location", "north");
    cb.Box="off";
    cb.TickDirection='out';
    cb.Label.String = "Magnitude";
    title(t,'LFP, MUA histogram, and first spike time relationship')
    cleanAxis(axs);
    linkaxes(axs,'x')
    saveFigure(f,fullfile(spk_file.folder,'Figures', ...
        'LFP - FST median relationship'),1,0)
    close all
end