function waveforms = getClusterWaveforms64Channels( bin_path, spike_subs)


cwf = zeros( 64, spikeWaveTime, numel( spkSubs ) );
frewind( fID );
for cspk = 1:numel(spkDists)
fseek( fID, 2*64*(spkDists(cspk)-spikeSamples), "cof" );
cwf(:,:,cspk) = fread( fID, [64, spikeWaveTime], "int16=>single" );
fseek( fID, -2*64*(spikeSamples+1), 'cof' );
end

end