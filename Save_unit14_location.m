load('Z:\Emilio\SuperiorColliculusExperiments\Roller\PoolFigures\eOPN3_ephMI.mat')
poolEphMI
lPSTH
poolEphMI
clearvars
load('Z:\Emilio\SuperiorColliculusExperiments\Roller\Batch18_ephys\MC\GADi43\240227_C+F_2200\ephys_E1\Unit14_localization.mat')

chanMap = load('Z:\Emilio\Cambridge_NeuroTech_F.mat');
int2cn = load('Z:\Emilio\IntanAmp2CambridgeNeurotech_Mapping.mat');
[shM, shB] = lineariz(chanMap.xcoords, 6, 1);
figure; plot(chanMap.xcoords, chanMap.xcoords.^[1,0] * [shM; shb])
chanMap.xcoords
figure; plot(chanMap.xcoords, chanMap.xcoords'.^[1,0] * [shM; shb])
figure; plot(chanMap.xcoords, chanMap.xcoords'.^[1,0] * [shM; shB])
figure; plot(chanMap.xcoords, chanMap.xcoords'.^[1,0] * [shM; shB], ".")
shMem = chanMap.xcoords'.^[1,0] * [shM;shB];
figure; plot(shMem, round(shMem))
figure; plot(shMem, round(shMem), '.')

fr
fs
whos fs
whos fr
save(fullfile(dataDir, 'GADi18_decode_sync.mat'), 'spike_times','velocity_GLM', 'fr', 'fs', 'gclID')
spike_times
velocity_GLM = velocity_GLM(:);
save(fullfile(dataDir, 'GADi18_decode_sync.mat'), 'spike_times','velocity_GLM', 'fr', 'fs', 'gclID')
save(fullfile(dataDir, 'Variables220823.mat'))
%-- 8/24/2022 14:32 --%
load('Z:\Emilio\SuperiorColliculusExperiments\Anaesthetised\EphysFigLaserAndWhiskers v2.mat')
sum(rclIdx)

gch = setdiff( chanMap, [50, 16, 17, 20] );
cch = sort( gch(mindist(1:13)) );

probe_geometry = [xcoords(gch), ycoords(gch)];
ptp = squeeze( range( cwf_norm(cch,:,:), 2 ) );

obj_func = @(idx,theta) ptp(:,idx) - (theta(1) ./ sqrt( sum( (probe_geometry - theta([2,3])).^2, 2 ) + theta(4).^2 ) );
loc_hat = zeros( size( ptp, 2 ), 4 );
parfor x = 1:size(ptp, 2)
    loc_hat(x,:) = fminsearch(@(w) sum( obj_func(x, w).^2 ), ...
        theta_init, optimset('Display', 'final', 'MaxIter', 4*2e4) );
end
loc_hat = zeros( size( ptp, 2 ), 4 ); cnvg_flag = false( size( ptp, 2 ), 1 );
parfor x = 1:size(ptp, 2)
    [loc_hat(x,:), ~, cnvg_flag(x)] = fminsearch(@(w) sum( obj_func(x, w).^2 ), ...
        theta_init, optimset('Display', 'none', 'MaxIter', 4*2e4) );
end
loc_hat
cnvg_flag
figure; [bcnt, xbedgs, ybedgs] = histcounts2( loc_hat(cnvg_flag,2), loc_hat(cnvg_flag,3), "NumBins", [128, 128], "XBinLimits", [-100, 1.2*max( xcoords )], "YBinLimits", [-100, 1.2*max( ycoords )] );
figure; imagesc( mean( [xbedgs(1:end-1);xbedgs(2:end)] ), mean( [ybedgs(1:end-1);ybedgs(2:end)] ), log10( bcnt +1 )', "AlphaData", bcnt )
figure; imagesc( mean( [xbedgs(1:end-1);xbedgs(2:end)] ), mean( [ybedgs(1:end-1);ybedgs(2:end)] ), log10( bcnt +1 )' )
axis xy
hold on; line( xcoords, ycoords, 'LineStyle', 'none', 'Marker', 'square', 'Color', 'k', 'MarkerFaceColor', 'k' )
figure; line( (xcoords(cch) + linspace(-20, 20, 74))', (ycoords(cch) + 20*( mean( cwf_norm(cch,:,cnvg_flag)./max( abs( cwf_norm(cch,:,cnvg_flag) ), [], 2 ), 3 ) ) )', 'color', 'k')
figure; line( (xcoords + linspace(-20, 20, 74))', (ycoords + 20*( mean( cwf_norm(:,:,cnvg_flag)./max( abs( cwf_norm(:,:,cnvg_flag) ), [], 2 ), 3 ) ) )', 'color', 'k')
axis square
cb = colorbar("Box", "off", "Location", "manual", "Color", "none");
cb = colorbar("Box", "off", "Location", "manual", "Color", "none", "TickDirection", "out", "AxisLocation", "out" );
cb = colorbar("Box", "off", "Location", "manual", "Color", 0.85*ones(1,3), "TickDirection", "out", "AxisLocation", "out" );
cb = colorbar("Box", "off", "Location", "south", "Color", 0.85*ones(1,3), "TickDirection", "out", "AxisLocation", "out" );
cb.Label.String = 'Log-Likelihood';
yticks( 2200 - yticks )
yticklabels( 2200 - yticks )
ylabel('Ventral \leftrightarrow Dorsal [mm]')
xlabel('Medial \leftrightarrow Lateral [mm]')
set( gca, "TickDir", "out" )
set( gca, "Box", "off", "Color", "none" )
xticklabels( 3600 - xticks )
xticklabels( 3600 + range(xcoords)/2 - xticks )
colormap( -inferno + 1 )
colormap( -jet + 1 )
colormap( -parula + 1 )
colormap( viko )
colormap( turbo )
colormap( -gray + 1 )
colormap( -bone + 1 )
colormap( parula )
colormap( -reds + 1 )
colormap( -greens + 1 )
colormap( -blues + 1 )
colormap( traffic )
colormap( -gray + 1 )
colormap( parula )
colormap( viko )
colormap( inferno )
colormap( -inferno + 1 )
colormap( inferno )
hold on; line( xcoords, ycoords, 'LineStyle', 'none', 'Marker', 'square', 'Color', 'k', 'MarkerFaceColor', 0.85*ones(1,3) )
colormap( turbo )
colormap( inferno )
colormap( turbo )
xticklabels( xticks - 3600 - range(xcoords)/2  )
xticklabels( xticks - 3600 + range(xcoords)/2  )
xticks
xticklabels( xticks + 3600 - range(xcoords)/2  )
title('Unit 14 location estimation')
xticklabels( xticks + 3600 - range(xcoords)/2  )
yticklabels( 2200 - yticks )
set( gca, "TickDir", "out" )
ylabel('Ventral \leftrightarrow Dorsal [mm]')
xlabel('Medial \leftrightarrow Lateral [mm]')
set( gca, "Box", "off", "Color", "none" )
saveFigure( figure(1), fullfile( "Z:\Emilio\SuperiorColliculusExperiments\Roller\Batch18_ephys\MC\GADi43\240227_C+F_2200\ephys_E1\Figures\Ephys VW-300.00-400.00 RW20.00-200.00 SW-300.00--120.00", "Unit14_location_estimation_map" ), true )
saveFigure( figure(2), fullfile( "Z:\Emilio\SuperiorColliculusExperiments\Roller\Batch18_ephys\MC\GADi43\240227_C+F_2200\ephys_E1\Figures\Ephys VW-300.00-400.00 RW20.00-200.00 SW-300.00--120.00", "Unit14_mean_waveforms_64ch" ), true )