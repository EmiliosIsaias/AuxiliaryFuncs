function cwf = getClusterWaveforms64Channels( bin_path, spike_subs)

cwf = [];
spikeWaveTime = 75;
spikeSamples = (spikeWaveTime - 1)/2;
tocol = @(x) x(:);
spkDists = [spike_subs(1); tocol( diff( spike_subs ) )];
fID = fopen( bin_path, "r" );
if fID >= 3
    cwf = zeros( 64, spikeWaveTime, numel( spike_subs ) );
    for cspk = 1:numel(spkDists)
        fseek( fID, 128*(spkDists(cspk)-spikeSamples), "cof" );
        cwf(:,:,cspk) = fread( fID, [64, spikeWaveTime], "int16=>single" );
        fseek( fID, -128*(spikeSamples+1), 'cof' );
    end
    fclose(fID);
end
end