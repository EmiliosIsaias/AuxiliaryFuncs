
fnOpts = {'UniformOutput', false};

home_path = '/gpfs/bwfor/home/hd/hd_hd/hd_bf154/';
repo_paths = cellfun(@(x) char( fullfile( home_path, x) ), ...
    {'NeuroNetzAnalysis', 'AuxiliaryFuncs', 'Scripts'},  fnOpts{:} );
addpath( repo_paths{:} )
%%
expandPath = @(x) fullfile( x.folder, x.name);
exp_path = getParentDir( pwd, 1);
m = 1e-3;

video_paths = dir( fullfile( exp_path, "Behavio*", "roller*.avi" ) );
vidObj = arrayfun(@(x) VideoReader( expandPath( x ) ), ...
    video_paths, fnOpts{:} );

readCSV = @(x) readtable(x, "Delimiter", ",");
fid_paths = dir( "FrameID*.csv" );

vidTx = arrayfun(@(x) readCSV( expandPath( x ) ), fid_paths, fnOpts{:} );
vidTx = cellfun(@(x) x.Var2 ./ 1e9, vidTx, fnOpts{:} ); % nanoseconds

dlc_paths = dir( fullfile( exp_path, "Behavio*", "roller*filtered.csv" ) );
fid_paths = dir( fullfile( exp_path, "Behavio*", "FrameID*.csv" ) );

tf_paths = dir( fullfile( exp_path, "Behavio*", "TriggerSignals*.bin") );
fsf_path = dir( fullfile( exp_path, "ephys*", "*_sampling_frequency.mat") );

if isempty( tf_paths )
    tf_paths = dir( fullfile( exp_path, "ephys*", "TriggerSignals*.bin") );
end

if isempty( fsf_path )
    fsf_path = dir( fullfile( exp_path, "Behavio*", "*_sampling_frequency.mat") );
end

fs_ephys = load( expandPath( fsf_path ), "fs" ); fs_ephys = fs_ephys.fs;
Ns_intan = [tf_paths.bytes]' ./ 4; % 2 signals x 2 bytes per sample.
Texp_ephys = Ns_intan ./ fs_ephys;

Texp_vid = cellfun(@(x) diff( x([1,end]) ), vidTx );

dlcTables = arrayfun(@(x) readDLCData(expandPath(x)), ...
    dlc_paths, fnOpts{:});

delta_tiv = Texp_ephys - Texp_vid;

lsrInt = cellfun(@(v,t) getLaserIntensitySignalFromVideo(v, t), ...
    vidObj(:), dlcTables(:), fnOpts{:} );

save( fullfile( video_paths(1).folder, "VideoLaserIntensity.mat" ), ...
    "lsrInt", "delta_tiv", "Texp_vid", "Texp_ephys")