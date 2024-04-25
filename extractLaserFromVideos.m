function [lsrInt, delta_tiv, Texp_vid, Texp_ephys] = ...
    extractLaserFromVideos( beh_path )

fnOpts = {'UniformOutput', false};
expandPath = @(x) fullfile( x.folder, x.name );
pathHere = @(x) fullfile( beh_path, x );
varsInFile = {'lsrInt', 'delta_tiv', 'Texp_vid', 'Texp_ephys'};

fprintf( 1, "Working directory: %s\n", beh_path )

if exist(pathHere( "VideoLaserIntensity.mat" ), "file")
    fprintf(1, "Laser signal file exists!\n")
    fprintf(1, "Loading... ")
    try
        load( pathHere( "VideoLaserIntensity.mat" ) , varsInFile{:} )
        fprintf(1, "Ready\n")
    catch
        fprintf(1, "Error while loading!\n")
    end
    return
else
    fprintf(1, "Computing laser signals from video(s)... \n")
end

exp_path = getParentDir( beh_path, 1);

m = 1e-3;

video_paths = dir( pathHere( "roller*.avi" ) );
vidObj = arrayfun(@(x) VideoReader( expandPath( x ) ), ...
    video_paths, fnOpts{:} );

readCSV = @(x) readtable(x, "Delimiter", ",");  
fid_paths = dir( pathHere( "FrameID*.csv" ) );

vidTx = arrayfun(@(x) readCSV( expandPath( x ) ), fid_paths, fnOpts{:} );
vidTx = cellfun(@(x) x.Var2 ./ 1e9, vidTx, fnOpts{:} ); % nanoseconds

dlc_paths = dir( pathHere( "roller*filtered.csv" ) );

tf_paths = dir( pathHere( "TriggerSignals*.bin") );
fsf_path = dir( fullfile( exp_path, "ephys*", "*_sampling_frequency.mat") );

if isempty( tf_paths )
    tf_paths = dir( fullfile( exp_path, "ephys*", "TriggerSignals*.bin") );
end

if isempty( fsf_path )
    fsf_path = dir( pathHere( "*_sampling_frequency.mat") );
end

fprintf(1, "Sampling frequency file: %s\n", fsf_path.name)

fs_ephys = load( expandPath( fsf_path ), "fs" ); fs_ephys = fs_ephys.fs;
Ns_intan = [tf_paths.bytes]' ./ 4; % 2 signals x 2 bytes per sample.
Texp_ephys = Ns_intan ./ fs_ephys;

Texp_vid = cellfun(@(x) diff( x([1,end]) ), vidTx );

dlcTables = arrayfun(@(x) readDLCData(expandPath(x)), ...
    dlc_paths, fnOpts{:});

delta_tiv = Texp_ephys - Texp_vid;

lsrInt = cellfun(@(v,t) getLaserIntensitySignalFromVideo(v, t), ...
    vidObj(:), dlcTables(:), fnOpts{:} );

save( pathHere( "VideoLaserIntensity.mat" ), varsInFile{:} )

end