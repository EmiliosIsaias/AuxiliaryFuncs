expandName = @(x) fullfile( x.folder, x.name );
roller_path = "Z:\Emilio\SuperiorColliculusExperiments\Roller";
figure_path = fullfile( roller_path, 'PoolFigures' );
owfFlag = false;
m = 1e-3; k = 1e3;
fnOpts = {'UniformOutput', false};
cellcat = @(x,d) cat( d, x{:} );
tocol = @(x) x(:);
%%
distPercent = 0.25;
pePaths = ["Batch12_ephys.e\MC\vGlut1\221206_C+F_2100",...
    "Batch11_ephys.MC\ChR2\Rb69\221128_F+C_2370", ...
    "Batch17_ephys.MC\eOPN3\Rb27\231103_C_2500", ...
    "Batch10_ephys\MC\WTg64\221111_C+PTX10microM_1900", ...
    "Batch10_ephys\MC\WTg64\221110_C_1900"];
vgePaths = ["Batch18_ephys\MC\GADi41\240213_F+C_2450", ...
    "Batch16_ephys\MC\GADi21\231023_C+F_2350", ...
    "Batch12_ephys.e\MC\vGlut5\221210_F+C_2450", ...
    "Batch11_ephys.MC\eOPN3\WT67\221202_C_2240", ...
    "Batch11_ephys.MC\ChR2\Rb71\221130_F+C_2100", ...
    "Batch11_ephys.MC\ChR2\Rb69\221129_C+F+Mus_2000", ...
    "Batch10_ephys\MC\WTg65\221115_E+PTX20microM_1900", ...
    "Batch8_ephys\MC\GADe55\220825_C+F_DV1750", ...
    "Batch7_ephys\MC\GADi52\220809_C+F_1900"];
aePaths = cat( 1, pePaths', vgePaths' );
Nexp = numel( aePaths );
rstVars2load = {'relativeSpkTmsStruct', 'configStructure'};
afVars2load = {'Conditions', 'fs'};
%%
vWin = [-300,400]*m;
clInfo = arrayfun(@(x) getClusterInfo( expandName( ...
    dir( fullfile( roller_path, x, '*\cluster_info.tsv' ) ) ) ), ...
    aePaths, fnOpts{:} );
%%
Nu = cellfun(@(x) sum( x.ActiveUnit ), clInfo );
Nuinit = cumsum( [1; Nu(1:end-1) ] );
Nuend = cumsum( Nu );
PSTHall_mu = zeros( 700, sum( Nu ) );
brAll = [];
PSTHall = cell( numel( aePaths ), 1 );
%%
for ce = 1:Nexp
    data_dir = fullfile( roller_path, aePaths(ce) );
    % condStruct = load( expandName( dir( fullfile( data_dir, "*\*analysis.mat" ) ) ), afVars2load{:} );
    rstPath = dir( fullfile( data_dir, "*\*RW20.00-200.00*RelSpkTms.mat" ) );
    brPath = dir( fullfile( data_dir, "*\BehaviourResult*.mat" ) );
    brVars2load = 'behRes';
    if isempty( brPath )
        brPath = dir( fullfile( data_dir, "*\Simple summary.mat" ) );
        brVars2load = 'summStruct';
    end
    brStruct = load( expandName( brPath ), brVars2load );
    behRes = brStruct.(brVars2load);
    ctrlSub = ismember( string( {behRes.ConditionName} ), 'Control Puff' );
    brAll = cat( 1, brAll, behRes(ctrlSub) );
    if ~isempty(rstPath) || numel( rstPath ) == 1
        rstCont = load( expandName( rstPath ), rstVars2load{:} );
    else
        fprintf(1, 'Either empty or more than 1 file found!\n');
        disp( {rstPath.name} )
        continue
    end
    rstStruct = rstCont.relativeSpkTmsStruct;
    confStruct = rstCont.configStructure;
    % Conditions = condStruct.Conditions;
    % fs = condStruct.fs;
    if any( confStruct.Viewing_window_s ~= vWin )
        confStruct.Viewing_window_s = vWin;
    end
    [PSTH, trial_tx, Na] = getPSTH_perU_perT( rstStruct, confStruct );
    a = Nuinit(ce); b = Nuend(ce);
    idx = a:b;
    PSTHall_mu(:,idx) = squeeze( mean( PSTH{1}, 1 ) );
    PSTHall(ce) = PSTH(1);
end
clearvars PSTH brStruct ctrlSub rstStruct confStruct brPath brVars2load ...
    rstCont a b idx;
ai_pt = arrayfun(@(s) getAIperTrial( s ), brAll, fnOpts{:} );
%%
[coeff, score, ~, ~, explained] = pca( zscore( ...
    PSTHall_mu(my_xor( trial_tx < [0,350]*m),:) ) );

y = pdist( coeff(:,1:3), 'euclidean' );
z = linkage( y, 'ward' );
rtm = cluster( z, 'cutoff', max(z(:,3))*distPercent, 'criterion', 'distance' );
Ncl = max( unique( rtm ) );
rl_mu = arrayfun(@(x) mean( zscore( PSTHall_mu(:,rtm==x) ), 2 ), ...
    unique(rtm), fnOpts{:} );
rl_mu = cat( 2, rl_mu{:} );

possCols = [2,3,5]; 
divCol = find( mod(Ncl,possCols) == 0, 1, "first" );
if isempty(divCol)
    divCol = ceil( sqrt( Ncl ) );
    Ncols = divCol;
else
    Ncols = possCols( divCol );
end
Nrows = ceil( Ncl/Ncols );
%%
respWins = [20,200; % Whole window
    20,50;          % sensory window
    50,200] * m;    % motor window in milliseconds
sponWins = -flip( respWins, 2 );
getBoolWindows = @(w) arrayfun(@(x) my_xor( trial_tx > w(x,:) ), ...
    1:size( w, 1 ), fnOpts{:} );

respFlags = cellcat( getBoolWindows( respWins ), 2 );
sponFlags = cellcat( getBoolWindows( sponWins ), 2 );
%%
% Response validation
ruFlags = false( sum( Nu ), 3 );
for cp = 1:Nexp
    auxP = false( Nu(cp), 3 );
    for cu = 1:Nu(cp)
        for cr = 1:size( respWins, 1 )
            d1 = squeeze( sum( PSTHall{cp}(:,sponFlags(:,cr),cu), 2 ) );
            d2 = squeeze( sum( PSTHall{cp}(:,respFlags(:,cr),cu), 2 ) );
            [~,auxP(cu,cr)] = poisson_means_etest( d1, d2, 0.05 );
        end
    end
    ruFlags(Nuinit(cp):Nuend(cp),:) = auxP;
end
clearvars auxP cp cu cr;

%%

[~, peak_resp] = max( rl_mu(respFlags(:,1), : ) - ...
    median( rl_mu(sponFlags(:,1),:) ), [], 1 );
elFlags = arrayfun( @(r) my_xor( k*respWins(1,1) + peak_resp(:) > ...
    respWins(r,:)*k ), [2,3], fnOpts{:} );
elFlags = cat( 2, elFlags{:} );
sponZSum = cellcat( arrayfun(@(r) sum( rl_mu(sponFlags(:,r),:) )', ...
    1:3, fnOpts{:} ), 2 );
respZSum = cellcat( arrayfun(@(r) sum( rl_mu(respFlags(:,r),:) )', ...
    1:3, fnOpts{:} ), 2 );

ctf = sqrt( (sponZSum(:,2:3) - respZSum(:,2:3)).^2 );
y2 = pdist( ctf, "cosine" );
z2 = linkage( y2, "complete" );
ctm = cluster( z2, 'cutoff', max(z2(:,3))*distPercent, ...
    'criterion', 'distance');
Ncl2 = max( unique( ctm ) );
%%
Nt_pe = cellfun('size', PSTHall, 1);
eelVals = zeros( sum( Nt_pe ), Ncl2 );
for cctm = 1:Ncl2
    init = 1;
    crtm = find( ctm == cctm );
    Nu_flag = any(rtm == crtm', 2);
    for cp = 1:Nexp
        eelVals(init:sum(Nt_pe(1:cp)), cctm) = ...
            mean( PSTHall{cp}(:, respFlags(:,1), ...
            Nu_flag(Nuinit(cp):Nuend(cp))), [2, 3] );
        init = init + Nt_pe(cp);
    end
end
expID = cellcat(arrayfun(@(x) x+ones(Nt_pe(x+1),1), 0:(numel(Nt_pe)-1), ...
    fnOpts{:} ), 1 );
%%
clrMap = ember( Nexp );
titls = ["Early", "Late", "Early only"];
rsPop = zeros( Nexp, Ncl2 );
for cctm = 1:Ncl2
    f = figure("Color", "w"); t = createtiles(f, 1, 1); ax = nexttile(t);
    axOpts = cleanAxis(ax); hold(ax, "on" );
    arrayfun(@(x) scatter(ax, zscore(eelVals(expID == x, cctm)), ...
        zscore( ai_pt{x} ), 32, clrMap(x,:), "filled", "o", ...
        "MarkerEdgeColor", clrMap(end-x+1,:), 'MarkerFaceAlpha', 0.7 ), ...
        1:Nexp ); title(ax, titls(cctm) )
    xlabel(ax, 'Lower \leftarrow SC activity \rightarrow Higher')
    ylabel(ax, 'Lower \leftarrow Startle response \rightarrow Higher')
    set( ax, 'TickDir', 'out' )
    lmObjs = arrayfun(@(x) fitlm( zscore(eelVals(expID == x, cctm)), ... 
        zscore( ai_pt{x} ), 'poly1' ), 1:Nexp, fnOpts{:} );
    rsPop(:,cctm) = cellfun(@(x) x.Rsquared.Ordinary, lmObjs );
    rfName = sprintf("%s activity vs startle", titls(cctm) );
    rfPath = fullfile( figure_path, rfName );
    if ~exist( rfPath + ".fig", "file" )
        saveFigure( f, rfPath, true, owfFlag )
    end
end
%%
f = figure("Color", "w"); t = createtiles(f,1,1); ax = nexttile(t); 
boxchart(ax, rsPop, "Notch", "on","BoxFaceColor","k", "MarkerStyle","none");
hold( ax, "on" )
swarmchart( ax, tocol( ones( size( rsPop, 1), 1) * (1:3) ), rsPop(:), ...
    'black','filled','o');
line( ax, (1:3)' * ones( 1, size( rsPop, 1) ), rsPop', 'Color', 'k', ...
    'LineWidth', 0.5 )
set( ax, 'TickDir', 'out' )
xticklabels( ax, titls )
ylabel(ax,'R²')
title(ax, 'R² for SC activity vs startle response')
pfName = "R² for SC activity vs startle response";
pfPath = fullfile( figure_path, pfName );
if ~exist( pfPath + ".fig", "file" )
    saveFigure( f, pfPath, true, owfFlag )
end
%%
f = figure('Color', 'w'); t = createtiles( f, Nrows, Ncols);
ax = gobjects( max( unique( rtm ) ), 1 );
clrMap = ember( Ncl );
for cm = unique( rtm )'
    ax(cm) = nexttile( t ); hold( ax(cm), 'on' )
    line( ax(cm), trial_tx*k, rl_mu(:,cm), 'Color', clrMap(cm,:), 'LineWidth', 2 )
    % line( ax(cm), k*trial_tx, movsum(rl_mu(:,cm), [10,0] ), ...
    %     'Color', clrMap(cm,:) )
    title( ax(cm), sprintf('Component %d, Members %d', cm, sum( rtm == cm ) ) )
    if ceil( cm / Ncols ) == 1
        disappearAxis( ax(cm), 'XAxis' );
    end
    if mod( cm, Ncols ) == 1
        ylabel(ax(cm), 'Z-score' )
    end
end
xlabel( ax, 'Time [ms]' )
linkaxes( ax, 'x' ); set( ax, 'TickDir', 'out' ); ytickangle(ax, 90)
xlim(ax(end), [-100, 350] )
phName = "PSTH PCA Classified";
phPath = fullfile( figure_path, phName );
if ~exist( phPath + ".fig", "file" )
    saveFigure( f, phPath, true, owfFlag )
end
%%
