function mean_delay = alignVideoWithEphys( lsrInt, trig, fs, beh_path )

expandPath = @(x) fullfile( x.folder, x.name);
fnOpts = {'UniformOutput', false};

% tf_paths = dir( fullfile( beh_path, "TriggerSignals*.bin") );
% fsf_path = dir( fullfile( exp_path, "ephys*", "*_sampling_frequency.mat") );
rs_path = dir( fullfile( beh_path, "RollerSpeed*.mat" ) );
load( expandPath( rs_path ), "fr")


% fIDs = arrayfun(@(x) fopen( expandPath( x ), "r" ), tf_paths );
% trig = arrayfun(@(x) fread( x, [2, inf], "uint16=>uint16" ), fIDs , ...
    % fnOpts{:} );
trig = cellfun(@(x) x', trig, fnOpts{:} );
% [~] = arrayfun(@(x) fclose( x ), fIDs );

lsrInt = cellfun(@(c) c - movmedian( c, round( 3*fr ) ), ...
    lsrInt, fnOpts{:} );
mean_delay = zeros( numel( lsrInt ), 1, "double" );
parfor cli = 1:numel(lsrInt)
    swObj = StepWaveform( trig{cli}(:,2), fs, 'verbose', false );
    testSubs = swObj.subTriggers;
    if numel(testSubs)
        lsrInt_loop = interp1( (0:length(lsrInt{cli})-1)' / fr, ...
            zscore( lsrInt{cli}(:) ), (0:length(trig{cli})-1)/fs, ...
            "pchip", "extrap");
        lsrInt_loop = iirCombFilter( lsrInt_loop, fs, 'Q', 17.5, 'W0', 9 );
        lsrInt_loop = iirCombFilter( lsrInt_loop, fs, 'Q', 35, 'W0', 50 );
        [r, lags] = xcorr( zscore( single( trig{cli}(:,2) ) ), lsrInt_loop(:), ...
            round( fs ) , "normalized" );
        [~, mean_delay_sub] = max( r );
        mean_delay(cli) = double(lags( mean_delay_sub ))/fs;
    end
end


end