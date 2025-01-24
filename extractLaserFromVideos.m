function [lsrInt, delta_tiv, Texp_vid, Texp_ephys, vidTx, trig, ...
    dlcTables, fs_ephys] = extractLaserFromVideos( beh_path )

fnOpts = {'UniformOutput', false};
expandPath = @(x) fullfile( x.folder, x.name );
pathHere = @(x) fullfile( beh_path, x );
varsInFile = {'lsrInt', 'delta_tiv', 'Texp_vid', 'Texp_ephys'};
out_path = pathHere( "VideoLaserIntensity_shuffle2.mat" );

fprintf( 1, "Working directory: %s\n", beh_path )

readCSV = @(x) readtable(x, "Delimiter", ",");
exp_path = getParentDir( beh_path, 1);

tf_paths = dir( pathHere( "TriggerSignals*.bin") );
fsf_path = dir( fullfile( exp_path, "ephys*", "*_sampling_frequency.mat") );

dlc_paths = dir( pathHere( "roller*shuffle2*filtered.csv" ) );

if isempty( tf_paths )
    tf_paths = dir( fullfile( exp_path, "ephys*", "TriggerSignals*.bin") );
end

if isempty( fsf_path )
    fsf_path = dir( pathHere( "*_sampling_frequency.mat") );
end

fprintf(1, "Sampling frequency file: %s\n", fsf_path.name)

fid_paths = dir( pathHere( "FrameID*.csv" ) );



fIDs = arrayfun(@(x) fopen( expandPath( x ), "r" ), tf_paths );
trig = arrayfun(@(x) fread( x, [2, inf], "uint16=>uint16" ), fIDs, fnOpts{:} );
arrayfun(@fclose, fIDs);

Ndlc = numel(dlc_paths);
dlcTables = cell( Ndlc, 1 );
vidTx = cell( Ndlc, 1 );
if Ndlc
    parfor cdlc = 1:Ndlc
        dlcTables{cdlc} = readDLCData( expandPath( dlc_paths(cdlc) ) );
        if ~isempty( fid_paths ) && numel( fid_paths ) == Ndlc
            vidTx{cdlc} = readCSV( expandPath( fid_paths(cdlc) ) ) ;
            vidTx{cdlc} = vidTx{cdlc}.Var2 * 1e-9; % nanoseconds
        end
    end
else
    fprintf(1, 'No DLC data found!\n')
    return
end
if isempty( vidTx )
    fprintf(1, 'No FrameID files in folder\n')
end

if ~isempty(fsf_path)
    fs_ephys = load( expandPath( fsf_path ), "fs" );
    fs_ephys = fs_ephys.fs;
end

if ~exist("fs_ephys", "var") || isempty(fs_ephys) || fs_ephys == 0
    fs_ephys = 3e4;
end

if exist( out_path, "file" )
    fprintf(1, "Laser signal file exists!\n")
    fprintf(1, "Loading... ")
    try
        load( out_path , varsInFile{:} )
        fprintf(1, "Ready\n")
    catch
        fprintf(1, "Error while loading!\n")
    end
    return
else
    fprintf(1, "Computing laser signals from video(s)... \n")
end

video_paths = dir( pathHere( "roller*.avi" ) );
hpcFlag = false;
if ~strcmpi( computer, 'PCWIN64' )
    vidObj = arrayfun(@(x) VideoReader( expandPath( x ) ), ...
        video_paths, fnOpts{:} );
else
    fprintf(1, 'Will not extract laser from video...\n')
    hpcFlag = true;
end

Ns_intan = [tf_paths.bytes]' ./ 4; % 2 signals x 2 bytes per sample.
Texp_ephys = Ns_intan ./ fs_ephys;
if ~isempty(vidTx) && all( cellfun(@(c) ~isempty( c ), vidTx ) )
    Texp_vid = cellfun(@(x) diff( x([1,end]) ), vidTx );
    delta_tiv = Texp_ephys - Texp_vid;
else
    fprintf(1, 'No FrameID files in this directory!\n')
    Texp_vid = 0; delta_tiv = 0;
end

testObj = cellfun(@(c) StepWaveform(c(2,:), fs_ephys, ...
    'verbose', false ), trig);
testSubs = arrayfun(@(x) x.subTriggers, testObj, fnOpts{:} );
laser_flag = ~cellfun(@isempty, testSubs);

lsrInt = cell( numel( video_paths ), 1 );
if ~hpcFlag
    for cvid = find( laser_flag(:)' )
        %parfor cvid = 1:numel( vidObj )
        lsrInt{cvid} = getLaserIntensitySignalFromVideo(vidObj{cvid}, ...
            dlcTables{cvid});
    end
end
save( out_path, varsInFile{:} )

end