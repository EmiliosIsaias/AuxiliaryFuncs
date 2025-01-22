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

%% Easy analysis
rs_whole_exp = zeros( Nexp, 2 );
for cexp = 1:Nexp
    x = zscore( mean( PSTHall{cexp}(:, my_xor( trial_tx < [20, 50]*m ), : ), [2,3] ) );
    y = zscore( mean( PSTHall{cexp}(:, my_xor( trial_tx < [50, 200]*m ), : ), [2,3] ) );
    z = zscore( ai_pt{cexp} )';
    f = figure("Color", "w" ); t = createtiles( f, 4, 2 );
    ax(1) = nexttile( t );
    aux_sens = fitlm( x, z, 'poly1' );
    plot( ax(1), aux_sens );
    xlabel( ax(1), '20-50 ms'); ylabel( ax(1), 'Startleness' )
    ax(2) = nexttile( t );
    aux_motr = fitlm( y, z, 'poly1' );
    plot( ax(2), aux_motr );
    xlabel( ax(2), '50-200 ms'); ylabel( ax(2), 'Startleness' )
    lgObjs = findobj(f,  'Type', 'Legend' );
    set( lgObjs, 'Box', 'off' )
    ax(3) = nexttile( t, [2,2] );
    aux_sm_mdl = fitlm( [x, y], z, "linear" );
    plot( ax(3), aux_sm_mdl );
    xlabel( ax(3), sprintf( '20-50 ms %.3f R²', aux_sens.Rsquared.Ordinary ) );
    ylabel( ax(3), sprintf( '50-200 ms %.3f R²', aux_motr.Rsquared.Ordinary ) );
    zlabel( ax(3), 'Startleness' );
    title( ax(3), sprintf("R² %.3f", aux_sm_mdl.Rsquared.Ordinary ) )
    title( t, sprintf( 'Experiment %d', cexp ) )
    cleanAxis( ax );
    grid( ax(3), "on" );
    rs_whole_exp(cexp, :) = [aux_sens.Rsquared.Ordinary, ...
        aux_motr.Rsquared.Ordinary];
end
%% Sliding window analysis for behaviour correlation
% Idea is to slide a time window per unit per trial for getting an R²
slid_win_length = 20*m; time_slide = 5*m;
time_init = -50*m; time_stop = 400*m;
Nrs = (time_stop - time_init - slid_win_length) / time_slide;
r_squared = cell( Nexp, 1 );
Nt = cellfun( "size", PSTHall, 1 );
parfor cexp = 1:Nexp
    r_squared{cexp} =  zeros( Nu(cexp), Nrs );
    for cu = 1:Nu(cexp)
        cw = time_init + [0, slid_win_length];
        aux_rs = zeros( 1, Nrs ); ci = 1;
        aux_trial = cellcat( arrayfun(@(t) conv( PSTHall{cexp}(t, :, cu ), ...
            gausswin( 5 ), "same" ), 1:Nt(cexp), fnOpts{:} ), 1 );
        while cw(2) <= time_stop
            act_mu = mean( aux_trial(:, my_xor( trial_tx < cw ) ) , 2 );
            if ( sum( act_mu == 0 ) / numel(act_mu ) ) < 0.4
                aux_mdl = fitlm( zscore( act_mu )', zscore( ai_pt{cexp} )', 'poly1' );
                aux_rs(ci) = aux_mdl.Rsquared.Ordinary;
            end
            cw = cw + time_slide; ci = ci + 1;
        end
        r_squared{cexp}(cu,:) = aux_rs;
    end
end
r_squared_cat = cat( 1, r_squared{:} );

%% Organising PSTHs by maximum R²
rsMdl = fit_poly( [1,Nrs], [time_init, time_stop - slid_win_length], 1 );
rsq_tx = ( (1:Nrs)' .^ [1,0] ) * rsMdl;

[mx_per_unit, ~] = max( r_squared_cat, [], 2 );
[~, rsq_mx_unit] = sort( mx_per_unit, "descend" );
f = figure("Color", "w" ); t = createtiles( f, 1, 2 );
ax(1) = nexttile(t); imagesc( ax(1), rsq_tx * k, [], ...
    r_squared_cat( rsq_mx_unit, :) );
cb(1) = colorbar( ax(1), "northoutside", "Box", "off", ...
    "TickDirection", "out" );
ax(2) = nexttile(t); imagesc( ax(2), trial_tx * k, [], ...
    ( PSTHall_mu(:, rsq_mx_unit) ./ max( PSTHall_mu(:, rsq_mx_unit) ) )' );
cleanAxis( ax );
colormap( inferno ); set( ax, 'tickdir', 'out' );
ylabel( ax(1), 'Units' ); disappearAxis(ax(2), 'YAxis' );
xlabel( ax,'Time [ms]' ); ytickangle(ax, 90 )
linkaxes( ax, 'xy' ); xlim( ax(1), k*[time_init, time_stop - slid_win_length])
%% Organising PSTHs by maximum R² and it's latency
rsMdl = fit_poly( [1,Nrs], [time_init, time_stop - slid_win_length], 1 );
rsq_tx = ( (1:Nrs)' .^ [1,0] ) * rsMdl;

[mx_per_unit, time_max] = max( r_squared_cat, [], 2 );
[~, rsq_mx_unit] = sort( mx_per_unit, "descend" );
[~, final_ord] = sortrows( [time_max, mx_per_unit], [1, 2], ...
    {'ascend', 'descend'} );
f = figure("Color", "w" ); t = createtiles( f, 1, 2 );
ax(1) = nexttile(t); imagesc( ax(1), rsq_tx * k, [], ...
    r_squared_cat( final_ord, :) );
cb(1) = colorbar( ax(1), "northoutside", "Box", "off", ...
    "TickDirection", "out" );
ax(2) = nexttile(t); imagesc( ax(2), trial_tx * k, [], ...
    zscore( PSTHall_mu(:, final_ord), 0, 1 )' );
cleanAxis( ax );
colormap( inferno ); set( ax, 'tickdir', 'out' );
ylabel( ax(1), 'Units' ); disappearAxis(ax(2), 'YAxis' );
xlabel( ax,'Time [ms]' ); ytickangle(ax, 90 )
linkaxes( ax, 'xy' ); xlim( ax(1), k*[time_init, time_stop - slid_win_length])

%% PCA for response type clustering
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
%% Grouping different neural response types
Nt_pe = cellfun('size', PSTHall, 1);
eelVals = zeros( sum( Nt_pe ), Ncl2 );
for cctm = 1:Ncl2
    init = 1;
    crtm = find( ctm == cctm );
    Nu_flag = any( rtm == crtm', 2);
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
rsPop = zeros( Nexp, Ncl );
for cctm = 1:Ncl
    f = figure("Color", "w"); t = createtiles(f, 1, 1); ax = nexttile(t);
    cleanAxis(ax); hold(ax, "on" );
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
        % saveFigure( f, rfPath, true, owfFlag )
    end
end
%% Not grouping and testing the effect of individual response types
Nt_pe = cellfun('size', PSTHall, 1);
eelVals = zeros( sum( Nt_pe ), Ncl );
for cctm = 1:Ncl
    Nu_flag = rtm == cctm;
    init = 1;
    for cp = 1:Nexp
        eelVals(init:sum(Nt_pe(1:cp)), cctm) = ...
            mean( PSTHall{cp}(:, respFlags(:,1), ...
            Nu_flag(Nuinit(cp):Nuend(cp))), [2, 3] );
        init = init + Nt_pe(cp);
    end
end
expID = cellcat(arrayfun(@(x) x+ones(Nt_pe(x+1),1), 0:(numel(Nt_pe)-1), ...
    fnOpts{:} ), 1 );
%% All classes in one figure
clrMap = ember( Ncl );
f = figure("Color", "w"); t = createtiles(f, 1, 1); ax = nexttile(t);
axOpts = cleanAxis(ax); hold(ax, "on" );

arrayfun(@(c) scatter(ax, zscore(eelVals(:, c)), ...
    zscore( cellcat( ai_pt, 2 ) ), 32, clrMap(c,:), "filled", "o", ...
    "MarkerEdgeColor", clrMap(end-c+1,:), 'MarkerFaceAlpha', 0.7 ), 1:Ncl );
titulo = "Response type vs startle response";
title(ax, titulo )
xlabel(ax, 'Lower \leftarrow SC activity \rightarrow Higher')
ylabel(ax, 'Lower \leftarrow Startle response \rightarrow Higher')
set( ax, 'TickDir', 'out' )
legend( ax, "Class " + string(1:Ncl)', axOpts{:}, 'Location', 'best', ...
    'NumColumns', Ncols);
rfName = titulo;
rfPath = fullfile( figure_path, rfName );
if ~exist( rfPath + ".fig", "file" )
    saveFigure( f, rfPath, true, owfFlag )
end

%% One class per figure
f = figure("Color", "w"); t = createtiles( f, Nrows, Ncols );
ax = gobjects( Nrows, Ncols );
for cc = 1:Ncl
    [c, r] = ind2sub( [Ncols, Nrows], cc );
    ax(r,c) = nexttile(t); cleanAxis(ax(r,c)); hold( ax(r,c), "on"	)
    arrayfun(@(e) scatter( ax(r,c), zscore( eelVals(expID==e,cc) ), ...
        zscore( ai_pt{e} ), 'o', 'MarkerFaceColor', 0.15*ones(1,3), ...
        'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.5 ), 1:Nexp )
    titulo = sprintf("Class %d (Neurons %d) vs startle response", cc, sum(rtm==cc) );
    title(ax(r,c), titulo )
    xlabel( ax(r,c), 'Lower \leftarrow SC activity \rightarrow Higher')
    ylabel( ax(r,c), 'Lower \leftarrow Amp Index \rightarrow Higher')
    yline( ax(r,c), 0, '--', 'Color', 0.5*ones(1,3) )
    xline( ax(r,c), 0, '--', 'Color', 0.5*ones(1,3) )
    set( ax(r,c), 'TickDir', 'out' )
    if r ~= Nrows
        disappearAxis( ax(r,c), 'XAxis' )
    end
    if c ~= 1
        disappearAxis( ax(r,c), 'YAxis' )
    end
end
title(t, sprintf('%d Experiments', Nexp ) )
linkaxes( ax, 'xy')

cfName = 'Classes activity vs startle response';
cfPath = fullfile( figure_path, cfName );
if ~exist( cfPath + ".fig", "file" )
    saveFigure( f, cfPath, true, owfFlag )
end
%% R² values

lmObjs = arrayfun(@(e) arrayfun(@(c) fitlm( zscore(eelVals(expID == e, c)), ...
    zscore(  ai_pt{e} ), 'poly1' ), 1:Ncl, fnOpts{:} ), 1:Nexp, fnOpts{:} );
rsPop = cellcat( cellfun(@(y) cellfun(@(x) x.Rsquared.Ordinary, y ), ...
    lmObjs, fnOpts{:} ), 1 );
%% Plot R² values
f = figure("Color", "w"); t = createtiles( f, 1, 1 );
ax = nexttile(t);
boxchart( ax, rsPop, 'Notch', 'on', 'BoxFaceColor', 'k', ...
    'MarkerStyle', 'none' );
xlabel(ax, 'Response types'); ylabel(ax, 'R²' );
hold( ax, "on" ); line( ax, ones( size( rsPop, 2 ), 1 ) * (1:14), rsPop', ...
    'Color', 0.15*ones(1,3), 'Linewidth', 1/3)
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
