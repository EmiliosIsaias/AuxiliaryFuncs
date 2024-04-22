
fnOpts = {'UniformOutput', false};

home_path = '/gpfs/bwfor/home/hd/hd_hd/hd_bf154/';
repo_paths = cellfun(@(x) char( fullfile( home_path, x) ), ...
    {'NeuroNetzAnalysis', 'AuxiliaryFuncs', 'Scripts'},  fnOpts{:} );
addpath( repo_paths{:} )
%%
expandPath = @(x) fullfile( x.folder, x.name);
beh_path = pwd;
fprintf( 1, "Working directory: %s\n", beh_path )
% beh_path = "Z:\Emilio\SuperiorColliculusExperiments\Roller\Batch18_ephys\iRNs\GADi39\240206_C+F_2399\Behaviour";
exp_path = getParentDir( beh_path, 1);

m = 1e-3;

video_paths = dir( fullfile( beh_path, "roller*.mp4" ) );
vidObj = arrayfun(@(x) VideoReader( expandPath( x ) ), ...
    video_paths, fnOpts{:} );

readCSV = @(x) readtable(x, "Delimiter", ",");  
fid_paths = dir( fullfile( beh_path, "FrameID*.csv" ) );

vidTx = arrayfun(@(x) readCSV( expandPath( x ) ), fid_paths, fnOpts{:} );
vidTx = cellfun(@(x) x.Var2 ./ 1e9, vidTx, fnOpts{:} ); % nanoseconds

dlc_paths = dir( fullfile( beh_path, "roller*filtered.csv" ) );
fid_paths = dir( fullfile( beh_path, "FrameID*.csv" ) );

tf_paths = dir( fullfile( beh_path, "TriggerSignals*.bin") );
fsf_path = dir( fullfile( exp_path, "ephys*", "*_sampling_frequency.mat") );

if isempty( tf_paths )
    tf_paths = dir( fullfile( exp_path, "ephys*", "TriggerSignals*.bin") );
end

if isempty( fsf_path )
    fsf_path = dir( fullfile( beh_path, "*_sampling_frequency.mat") );
end

fprintf(1, "Sampling frequency file: %s\n", fsf_path.name)

fs_ephys = load( expandPath( fsf_path ), "fs" ); fs_ephys = fs_ephys.fs;
Ns_intan = [tf_paths.bytes]' ./ 4; % 2 signals x 2 bytes per sample.
Texp_ephys = Ns_intan ./ fs_ephys;

Texp_vid = cellfun(@(x) diff( x([1,end]) ), vidTx );

dlcTables = arrayfun(@(x) readDLCData(expandPath(x)), ...
    dlc_paths, fnOpts{:});

delta_tiv = Texp_ephys - Texp_vid;
%% 
lsrInt = cellfun(@(v,t) getLaserIntensitySignalFromVideo(v, t), ...
    vidObj(:), dlcTables(:), fnOpts{:} );

save( fullfile( beh_path, "VideoLaserIntensity.mat" ), ...
    "lsrInt", "delta_tiv", "Texp_vid", "Texp_ephys")