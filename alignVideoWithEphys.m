expandPath = @(x) fullfile( x.folder, x.name);

load('VideoLaserIntensity.mat')

tf_paths = dir( fullfile( beh_path, "TriggerSignals*.bin") );
fsf_path = dir( fullfile( exp_path, "ephys*", ...
    "*_sampling_frequency.mat") );
rs_path = dir( fullfile( beh_path, "RollerSpeed*.mat" ) );
load( expandPath( rs_path ), "fr")


fIDs = arrayfun(@(x) fopen( expandPath( x ), "r" ), tf_paths );
trig = arrayfun(@(x) fread( x, [2, inf], "uint16=>uint16" ), fIDs , ...
    fnOpts{:} );
trig = cellfun(@(x) x', trig, fnOpts{:} );
[~] = arrayfun(@(x) fclose( x ), fIDs );

lsrInt_aux = lsrInt{1} - movmedian( lsrInt{1}, round( 3*fr ) );
lsrInt_aux = interp1( (0:length(lsrInt{1})-1)/fr, zscore( lsrInt_aux(:) ), ...
    (0:length(trig{1})-1)/fs, "pchip", "extrap");

[r, lags] = xcorr( zscore( single( trig{1}(:,2) ) ), lsrInt_aux(:), ...
    round( fs ) , "normalized" );