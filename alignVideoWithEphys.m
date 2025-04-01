function [mean_delay, lsrInt, fr] = alignVideoWithEphys( lsrInt, trig, fs, beh_path )

expandPath = @(x) fullfile( x.folder, x.name);
fnOpts = {'UniformOutput', false};


% tf_paths = dir( fullfile( beh_path, "TriggerSignals*.bin") );
% fsf_path = dir( fullfile( exp_path, "ephys*", "*_sampling_frequency.mat") );
rs_path = dir( fullfile( beh_path, "RollerSpeed*.mat" ) );
if isempty(rs_path)
    readCSV = @(x) readtable(x, "Delimiter", ",");
    flfa = @(x) fullfile(x.folder, x.name);
    search4This = @(x) dir(fullfile(beh_path, x));
    fFiles = search4This("FrameID*.csv");

    vidTx = arrayfun(@(x) readCSV(flfa(x)), fFiles, fnOpts{:});
    vidTx = cellfun(@(x) x.Var2/1e9, vidTx, fnOpts{:}); % nanoseconds
    estFr = cellfun(@(x) median( 1 ./ diff(x) ), vidTx);
    fr = mean( estFr );
    % save( expandPath(rs_path), "fr" )
else
    load( expandPath( rs_path ), "fr" )
end


% fIDs = arrayfun(@(x) fopen( expandPath( x ), "r" ), tf_paths );
% trig = arrayfun(@(x) fread( x, [2, inf], "uint16=>uint16" ), fIDs , ...
% fnOpts{:} );
trig = cellfun(@(x) x', trig, fnOpts{:} );
% [~] = arrayfun(@(x) fclose( x ), fIDs );
lsrInt_bp = lsrInt;
lsrInt = cellfun(@(c) c - movmedian( c, round( 3*fr ) ), ...
    lsrInt, fnOpts{:} );
mean_delay = zeros( numel( lsrInt ), 1, "double" );
laser_flag = ~cellfun(@isempty, lsrInt);
if any(laser_flag)

    lsrInt = cellfun(@(c) iirCombFilter( c, fr, 'Q', 17.5, 'W0', 9, ...
        'verbose', false ), lsrInt(laser_flag), fnOpts{:} );
    lsrInt = cellfun(@(c) iirCombFilter( c, fr, 'Q', 35, 'W0', 50, ...
        'verbose', false ), lsrInt, fnOpts{:} );
    clii = 1;
    for cli = find(laser_flag')
        swObj = StepWaveform( trig{cli}(:,2), fs, 'verbose', false );
        testSubs = swObj.subTriggers;
        if numel(testSubs)
            lsrInt_loop = interp1( (0:length(lsrInt{clii})-1)' / fr, ...
                zscore( lsrInt{clii}(:) ), (0:length(trig{cli})-1)/fs, ...
                "pchip", "extrap");
            [r, lags] = xcorr( zscore( single( trig{cli}(:,2) ) ), lsrInt_loop(:), ...
                round( fs * 0.3 )  , "normalized" );
            [~, mean_delay_sub] = max( r );
            mean_delay(cli) = double(lags( mean_delay_sub ))/fs;
        end
        clii = clii + 1;
    end
end

end