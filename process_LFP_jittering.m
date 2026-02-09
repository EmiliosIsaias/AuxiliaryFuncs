fnOpts = {'UniformOutput', false};
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
for ce = experiments_lfp
    trig_file = dir(fullfile(fullpath(ce),'**','*analysis.mat'));
    if numel(trig_file) > 1
        [~,bigger_file] = max([trig_file.bytes]);
        trig_file = trig_file(bigger_file);
    end
    fs_file = dir(fullfile(fullpath(ce),'**', '*_sampling_frequency.mat'));
    if numel(fs_file)>1
        fprintf(1,'Don''t know what to do yet')
    end
    load(fullpath(fs_file),'fs')
    load(fullpath(trig_file), "Conditions", "Triggers");
    lfp = Triggers.LFP; clearvars Triggers;
    lfp_ds = resample(lfp,(((1:length(lfp))-0.5)/fs),fs_low);
    tx = ((1:length(lfp_ds))-0.5)/fs_low;
    lfp = filtfilt(b,a,lfp_ds);
    temp = zscore(lfp);
    brain_state = fitSpline(tx, temp, 1, 0.65, 0.75, 1);
    nf_flag = brain_state <= 0;
    sObj = StepWaveform(nf_flag,fs_low);
    sObj.MinIEI = 0.5;
    nfs_triggers = sObj.Triggers/fs_low;
    % nf_signal = false(1,length(lfp_ds));
    % for cr = 1:size(nfs_triggers,1)
    %     nf_signal = nf_signal | ...
    %         (tx>=nfs_triggers(cr,1) & tx<=nfs_triggers(cr,2));
    % end
    wscc = contains({Conditions.name}, {'whisker','piezo'}, "IgnoreCase", true) & ...
        contains({Conditions.name}, 'Control');
    ws_triggers = Conditions(wscc).Triggers(:,1)/fs;
    on_stim_flag = any(ws_triggers' > nfs_triggers(:,1) &...
        ws_triggers' < nfs_triggers(:,2), 1);
    delta_t1 = ws_triggers(:,1)' - nfs_triggers(:,1);
    delta_t2 = ws_triggers(:,1)' - nfs_triggers(:,2);
end