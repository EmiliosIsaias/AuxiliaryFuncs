expandName = @(x) fullfile( x.folder, x.name );
roller_path = fullfile( "Z:\Emilio\SuperiorColliculusExperiments\Roller" );
% Z:\Emilio\SuperiorColliculusExperiments\Roller\Batch1_ephys\BC\GADi48\210708_C
data_dir = fullfile( roller_path, 'Batch1_ephys\BC\GADi48\210708_C' );
owfFlag = false;
m = 1e-3; k = 1e3; 
fnOpts = {'UniformOutput', false};
%%
load( expandName( dir( fullfile( data_dir, "*\*analysis.mat" ) ) ), ...
    "Conditions", "fs" )
rstPath = dir( fullfile( data_dir, "*\*RW20.00-200.00*RelSpkTms.mat" ) );
if ~isempty(rstPath) || numel( rstPath ) == 1
    load( expandName( rstPath ), 'relativeSpkTmsStruct', 'configStructure' )
else
    fprintf(1, 'Either empty or more than 1 file found!\n');
    disp( {rstPath.name} )
    return
end
figure_path = expandName( dir( fullfile( data_dir, "*\Figur*" ) ) );
%% PSTH and PCA
[PSTH, trial_tx, Na] = getPSTH_perU_perT( relativeSpkTmsStruct, configStructure );
Nu = size( PSTH{1}, 3 );
PSTHt = squeeze( mean( PSTH{1}, 1 ) );
f = figure("Color", "w"); t = createtiles( f, 4, 1 ); ax(1) = nexttile(t, [3, 1]);
imagesc( ax(1), k*trial_tx, 1:Nu, zscore( PSTHt )' ); colormap( inferno );
hold ( ax(1), 'on' ); 
xline( ax(1), 0, 'Color', 0.85*ones(1,3), 'LineStyle', '--');
patch( ax(1), [20,50,50,20], [0, 0, Nu, Nu]+0.5, 1, 'FaceAlpha', 0.1, ...
    'edgecolor', 'none', 'FaceColor', [0,0.51,1] )
patch( ax(1), [50,200,200,50], [0, 0, Nu, Nu]+0.5, 1, 'FaceAlpha', 0.1, ...
    'edgecolor', 'none', 'FaceColor', [1,0,0.51] )

xlabel(ax(1), 'Time [ms]'); disappearAxis( ax(1), 'YAxis' ); 
[coeff, score, ~, ~, explained] = pca( zscore( ...
    PSTHt(my_xor( trial_tx < [0,350]*m),:) ) );
ax(2) = nexttile(t); line( ax(2), 1:numel( explained), cumsum( explained ) )
ylim( ax(2), [0,100] ); yline( ax(2), sum( explained(1:3) ), 'k' )
set( ax, 'TickDir', 'out' ); arrayfun(@(x) cleanAxis(x), ax, fnOpts{:} );

distPercentage = 0.25;
%% Cluster based on PCA. Does it look alright?
y = pdist( coeff(:,1:3), 'euclidean' );
z = linkage( y, 'ward' );
figure; dendrogram( z, Nu, "ColorThreshold", distPercentage*max(z(:,3)) )
%% Hierarchical clustering
rtm = cluster( z, 'cutoff', max( z(:,3) )*distPercentage, ...
    'criterion', 'distance' );
possCols = [2,3,5];
Ncl = max( unique( rtm ) );
divCol = find( mod( Ncl, possCols ) == 0, 1, "first" );
if ~numel(divCol)
    divCol = ceil( sqrt( Ncl ) );
    Ncols = divCol;
else
    Ncols = possCols( divCol );
end
Nrows = ceil( Ncl / Ncols );
figure; silhouette( coeff(:,1:3), rtm )
%% PCA Components
f = figure("Color", "w"); t = createtiles( f, 1, 1 );
clrMap = inferno( Ncl );
ax = nexttile( t ); hold(ax, 'on'); axOpts = cleanAxis( ax );
arrayfun(@(x) line( ax, coeff(rtm==x,1), coeff(rtm==x,2), coeff(rtm==x,3), ...
    'Marker', 'o', 'MarkerFaceColor', clrMap(x,:), 'MarkerEdgeColor', ...
    clrMap(end-x+1,:), 'LineStyle', 'none' ), unique( rtm ) )
ytickangle( ax, 90 ); ax.View = [75, 25];
legend( 'show', axOpts{:}, 'Location', 'best', 'NumColumns', Ncols);
xlabel( ax, 'PC1' ); ylabel( ax, 'PC2' ); zlabel( ax, 'PC3' );
title(ax, '3 PC for response type classification' )
pcaFigName = 'PCA components classified';
pcaFigPath = fullfile( figure_path, pcaFigName );
if ~exist( pcaFigPath+".fig", "file" )
    saveFigure( f, pcaFigPath, true, owfFlag );
end
%% Plot mean z-score PSTH for response type membership
f = figure('Color', 'w'); t = createtiles( f, Nrows, Ncols);
ax = gobjects( max( unique( rtm ) ), 1 );
rl_mu = zeros( size( PSTHt, 1 ), Ncl );
for cm = unique( rtm )'
    ax(cm) = nexttile( t ); 
    rl_mu(:,cm) = mean( zscore( PSTHt(:,rtm==cm) ), 2 );
    line( ax(cm), trial_tx*1e3, rl_mu(:,cm), 'Color', 'k' )
    hold( ax(cm), "on" )
    arrayfun(@(x) patch( ax(cm), 1e3*trial_tx([1:end,end:-1:1]), ...
        zscore( PSTHt([1:end,end:-1:1],x) ), 'k', 'EdgeAlpha', 0.05 ), ...
        find( rtm == cm ) )
    title( ax(cm), sprintf( 'Component %d, Members# %d', cm, sum( rtm == cm ) ) )
    cleanAxis( ax(cm) );
    if ceil( cm / Ncols ) ~= Nrows
        disappearAxis( ax(cm), 'XAxis' );
    end
    if mod( cm, Ncols ) == 1
        ylabel( ax(cm), 'Z-score' )
    end
end
xlabel( ax, 'Time [ms]')
ytickangle( ax, 90 )
set( ax, 'TickDir', 'out' )
linkaxes( ax, 'x' )
xlim(ax(cm), [15,250] )
arrayfun(@(x) xline( x, [20,50,200], 'k--' ), ax )
figName = "PSTH PCA classification";
figPath = fullfile( figure_path, figName );
if ~exist(figPath+".fig","file") || owfFlag
    saveFigure( f, figPath, true)
end
%% Assigning each component to early only, early, and late classes.
respWins = [20,200; % Whole window
    20,50;          % sensory window
    50,200] * m;    % motor window in milliseconds
sponWins = -flip( respWins, 2 );
getBoolWindows = @(w) arrayfun(@(x) my_xor( trial_tx > w(x,:) ), ...
    1:size( w, 1 ), fnOpts{:} );
cellcat = @(d,x) cat( d, x{:} );

respFlags = cellcat( 2, getBoolWindows( respWins ) );
sponFlags = cellcat( 2, getBoolWindows( sponWins ) );