expandName = @(x) fullfile( x.folder, x.name );
roller_path = "Z:\Emilio\SuperiorColliculusExperiments\Roller";
figure_path = fullfile( roller_path, 'PoolFigures' );
owfFlag = false;
m = 1e-3; k = 1e3;
fnOpts = {'UniformOutput', false};
cellcat = @(x,d) cat( d, x{:} );
tocol = @(x) x(:);
getMI = @(x,d) diff(x, 1, d) ./ ...
    sum(x, d).*(sum(x,d)>0) + (1.*(sum(x,d)==0 | sum(x,d)< 1e-12));
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
    "Batch10_ephys\MC\WTg65\221114_C_2370", ...
    "Batch10_ephys\MC\WTg65\221115_E+PTX20microM_1900", ...
    "Batch8_ephys\MC\GADe55\220825_C+F_DV1750", ...
    "Batch7_ephys\MC\GADi52\220809_C+F_1900"];
aePaths = cat( 1, pePaths', vgePaths' );
Nexp = numel( aePaths );
rstVars2load = {'relativeSpkTmsStruct', 'configStructure'};
afVars2load = {'Conditions', 'fs'};
%%
vWin = [-300, 400]*m;
clInfo = arrayfun(@(x) getClusterInfo( expandName( ...
    dir( fullfile( roller_path, x, 'ephys*\cluster_info.tsv' ) ) ) ), ...
    aePaths, fnOpts{:} );
%% Pooling ephys and behaviour
% Initialising
Nu = cellfun(@(x) sum( x.ActiveUnit ), clInfo );
Nuinit = cumsum( [1; Nu(1:end-1) ] );
Nuend = cumsum( Nu );
PSTHall_mu = zeros( 700, sum( Nu ) );
brAll = [];
PSTHall = cell( numel( aePaths ), 1 );
uSig = PSTHall; uID = PSTHall;
uMod = PSTHall;
uMI = PSTHall;
for ce = 1:Nexp
    data_dir = fullfile( roller_path, aePaths(ce) );
    % condStruct = load( expandName( dir( fullfile( data_dir, "*\*analysis.mat" ) ) ), afVars2load{:} );
    rstPath = dir( fullfile( data_dir, "*", "*RW20.00-200.00*RelSpkTms.mat" ) );
    brPath = dir( fullfile( data_dir, "*","BehaviourResult*.mat" ) );
    mfPath = dir( fullfile( data_dir, "ephys*", "Results", "Res VW* RW20.00-200.00 ms SW*PuffAll.mat") );
    brVars2load = 'behRes';
    if isempty( brPath )
        brPath = dir( fullfile( data_dir, "*\Simple summary.mat" ) );
        brVars2load = 'summStruct';
    end
    brStruct = load( expandName( brPath ), brVars2load );
    behRes = brStruct.(brVars2load);
    ctrlSub = ismember( string( {behRes.ConditionName} ), 'Control Puff' );
    brAll = cat( 1, brAll, behRes(ctrlSub) );
    if ~isscalar( rstPath )
        nms = {rstPath.name}';
        % Response, spontaneous, view
        params = cellfun(@(x) str2double( x ), ...
            cellcat( cellfun(@(y) regexp( y, '\d+\.\d+', 'match' ), ...
            nms, fnOpts{:} ), 1 ) );
        params = params .* -[-1,-1,1,1,1,-1];
        min_subs = [];
        sel = 1;
        if numel( nms ) > 1
            params_diff = diff( params, 1, 1 );
            if ~nnz( params_diff(3:end) )
                % Equivalent results
                %ce = find( ce, 1, "first" );
                params = params(1,:);
                nms = nms(1);
            else
                % Spontaneous window messed up. Selecting the furthest and
                % equal to response window's length.
                % ce = find( ce );
                rWins = diff( params(:,1:2), 1, 2 );
                sWins = diff( params(:,3:4), 1, 2 );
                eq_flag = sWins == max( rWins );
                if all(eq_flag)
                    [~, min_subs] = min( params(:,3:4), [], 1 );
                    % ce = ce(min_subs(1)); %#ok<*FXSET>
                    params = params(min_subs(1),:);
                    nms = nms(min_subs(1));
                    sel = min_subs(1);
                else
                    % ce = ce(eq_flag);
                    params = params(eq_flag);
                    nms = nms(eq_flag);
                    sel = nms(eq_flag);
                end
            end
        end
        rstPath = rstPath(sel);
    end
    rstCont = load( expandName( rstPath ), rstVars2load{:} );
    if ~isscalar(mfPath)
        mfNms = {mfPath.name}';
        mfSel = contains( mfNms, sprintf( 'VW%.2f-%.2f', params(:,5:6) ) ) & ...
            contains( mfNms, sprintf( 'RW%.2f-%.2f', params(:,1:2) ) ) & ...
            contains( mfNms, sprintf( 'SW%.2f-%.2f', params(:,3:4) ) );
        if nnz(mfSel)==1
            mfPath = mfPath(mfSel);
        else
            fprintf(1, 'Big problems!\n')
        end
    end
    mftype = 1;
    mfVars2load = {'Results', 'gclID', 'Counts'};
    if isempty( mfPath )
        fprintf( 1, 'Res file not found!\n')
        mfPath = dir( fullfile( data_dir, "ephys*", "Results", "Map*.mat") );
        mfVars2load = {'keyCell', 'resMap'};
        if isempty( mfPath )
            fprintf( 1, 'Unable to load unit response information!\n')
            disp( mfPath )
            continue
        end
        mftype = 2;
    end
    mfStruct = load( expandName( mfPath ), mfVars2load{:} );
    if mftype == 1
        Results = mfStruct.Results; gclID = mfStruct.gclID;
        Counts = mfStruct.Counts(1,:); Counts = mean( cat( 3, Counts{:} ), 2 );
        configFlag = cellfun(@(x) ~isempty(x), regexp( {Results.Combination}, ...
            '1\s1\ssignrank', 'ignorecase' ) );
        uSig{ce} = Results(configFlag).Activity(1).Pvalues;
        uMod{ce} = Results(configFlag).Activity(1).Direction;
        uMI{ce} = getMI( Counts, 3 );
        uID{ce} = gclID;
    else
        resMap = mfStruct.resMap; keyCell = mfStruct.keyCell;
        configFlag = contains( keyCell(:,1), 'RW20.00-200.00' ) & ...
            contains( keyCell(:,4), 'Control Puff' );
        if sum( configFlag ) ~= 1
            fprintf( 1, 'Cannot process this keycell... \n')
            disp( keyCell )
        end
        uSig{ce} = resMap( keyCell{configFlag,:} );
        uMod{ce} = zeros( Nu(ce), 1 ); uMI{ce} = uMod{ce};
        uID{ce} = clInfo{ce}{clInfo{ce}.ActiveUnit==1, "cluster_id" };
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
            gausswin( 5 ), "same" ), 1:Nt(cexp), 'UniformOutput', false ), 1 );
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
%%
r_max_all = max( r_squared_cat, [], 2 );
[~, rsq_mx_unit] = sort( r_max_all, "descend" );
f = figure("Color", "w" ); t = createtiles( f, 1, 7 );
ax(1) = nexttile(t, [1,3]); imagesc( ax(1), rsq_tx * k, [], ...
    r_squared_cat( rsq_mx_unit, :) );
cb = colorbar( ax(1), "northoutside", "Box", "off", ...
    "TickDirection", "out" );
cb.Label.String = "R²";
ax(2) = nexttile(t,[1,3]); imagesc( ax(2), trial_tx * k, [], ...
    ( PSTHall_mu(:, rsq_mx_unit) ./ max( PSTHall_mu(:, rsq_mx_unit) ) )' );
title(ax(2), 'Single unit PSTHs ordered by descending R²')
mgn_max_all = max( PSTHall_mu( my_xor( trial_tx < [20,200]*m ), rsq_mx_unit ) );
ax(3) = nexttile( t ); line( ax(3), mgn_max_all, 1:sum(Nu), 'Color', 'k' );
set( get( ax(3), 'YAxis' ), 'Direction', 'reverse' )
cleanAxis( ax );
colormap( inferno ); set( ax, 'tickdir', 'out' );
ylabel( ax(1), 'Units' ); arrayfun(@(x) disappearAxis(x, 'YAxis' ), ax(2:end) );
xlabel( ax(1:2),'Time [ms]' ); ytickangle(ax, 90 )
xlabel( ax(3), 'Spikes/bin')
linkaxes( ax(1:2), 'xy' ); linkaxes( ax([1,3]), 'y' )
xlim( ax(1), k*[time_init, time_stop - slid_win_length])
fz_sub = find( r_max_all(rsq_mx_unit) == 0, 1, "first" );
arrayfun(@(x) yline( x, fz_sub, 'Color', 0.85*ones(1,3) ), ax(1:2) )
yline( ax(3), fz_sub, 'Color', [0.0314, 0.1882, 0.4196] )
ylim( ax, [0, sum(Nu)] + [1,-1]/2 )
%% Time resolved boxplots for all experiments
for td = 30
    slid_win_length = td*m; time_slide = td*m;
    time_init = -160*m; time_stop = 400*m;
    Nrs = floor( (time_stop - time_init - slid_win_length) / time_slide );
    r_squared_pexp = zeros( Nexp, Nrs );
    for cexp = 1:Nexp
        cw = time_init + [0, slid_win_length];
        for ci = 1:Nrs
            act_mu = mean( PSTHall{cexp}(:, my_xor( trial_tx < cw ),: ) , [2,3] );
            aux_mdl = fitlm( zscore( act_mu )', zscore( ai_pt{cexp} )', 'poly1' );
            r_squared_pexp(cexp,ci) = aux_mdl.Rsquared.Ordinary;
            cw = cw + time_slide;
        end
    end
    aux_mdl = fit_poly( [1, Nrs], [time_init, time_stop - slid_win_length], 1);
    b_tx = ( ( 1:Nrs )'.^[1,0] ) * aux_mdl;
    % r_squared_cat = cat( 1, r_squared_pexp{:} );
    %% Plotting results
    bxOpts = {'JitterOutliers', 'on', 'MarkerStyle', '.', 'MarkerColor', 'k',...
        'BoxFaceColor', 'k', 'BoxWidth', k*(time_slide)/2, ...
        'Notch', 'off' };

    f = figure('Color', 'w'); t = createtiles(f, 1, 1); ax = nexttile( t );
    boxchart(ax, tocol( ones(Nexp,1)*b_tx' * k), ...
        tocol(r_squared_pexp), bxOpts{:} )
    hold( ax, 'on');
    line(ax, k*b_tx, median( r_squared_pexp, 1 ), 'Color', 'k', 'LineWidth', 2 )
    xline( ax, [0,50,200], 'r--')
    xlabel( ax, 'Time [ms]' )
    xlim( ax, k*(b_tx([1,end]) + [-1;1]*slid_win_length/2) )
    set( ax, 'TickDir', 'out' )
    ytickangle( ax, 90 )
    ylabel( ax, 'R² per window' )
    cleanAxis( ax );
    title( ax, sprintf('Time-resolved_{%d ms} R² population (per experiment)', td ) )
    [p, tbl, stats] = kruskalwallis( r_squared_pexp, string( b_tx(:) * k ) );
    figure; tbl2 = multcompare( stats );
    nnz(tbl2(:,end) < 0.05 )
end
%% Population impact on regression
r_max = cellfun(@(x) max( x, [], 2 ), r_squared, fnOpts{:} );
r_max_all = cat(1, r_max{:} );
time_flags = my_xor( trial_tx < [20,200]*m );
pop_ax = 0:0.1:0.9;
%pop_ax = 1:-0.1:0.1;
r_squared_ppop = zeros( numel( pop_ax ), Nexp );
for cexp = 1:Nexp
    for cpop = pop_ax
        act_mu = mean( zscore( PSTHall{cexp}(:, time_flags, ...
            r_max{cexp} >= ( max( r_max{cexp} ) * cpop ) ) ), [2,3] );
        % act_mu = mean( PSTHall{cexp}(:, time_flags, ...
        %     r_max{cexp} <= ( max( r_max{cexp} ) * cpop ) ), [2,3] );
         % act_mu = mean( PSTHall{cexp}(:, time_flags, ...
         %    r_max{cexp} >= ( max( r_max{cexp} ) * cpop ) ), [2,3] );
        % PSTHtest = PSTHall{cexp}./max(PSTHall{cexp}, [], [2,3] );
        % act_mu = mean( zscore( PSTHtest(:, time_flags, ...
        %     r_max{cexp} >= ( max( r_max{cexp} ) * cpop ) ) ), [2,3] );
        aux_mdl = fitlm( zscore( act_mu )', zscore( ai_pt{cexp} )', 'poly1' );
        r_squared_ppop(round(10*cpop+1),cexp) = ...
            aux_mdl.Rsquared.Ordinary;
        % r_squared_ppop(round(-10*cpop+11),cexp) = ...
        %     aux_mdl.Rsquared.Ordinary;
    end
end
% Plotting results
f = figure("Color","w"); t = createtiles( f, 1, 1);
ax = nexttile( t );
boxchart( ax, r_squared_ppop', "BoxFaceColor", "k", "MarkerStyle", ".", "Notch", "off" )
set( ax, 'TickDir', 'out' ); cleanAxis( ax ); ytickangle( ax, 90 );
ylabel( ax, 'R² by considered population')
xlabel( ax, 'All cells \leftrightarrow Lower R² cells [%]' )
xticklabels( ax, 100*(pop_ax) )
title( ax, 'R² relationship to the considered population')
%% R² histogram
f = figure("Color", "w"); t = createtiles( f, 1, 1 );
ax = nexttile( t ); histogram( r_max_all, "BinLimits", [0, 0.5], ...
    "BinWidth", 0.01, "Normalization", "probability", "EdgeColor", "none", ...
    "FaceAlpha", 1/3 )
set( get( ax, 'YAxis'), 'Scale', 'log' ); xlim( ax, [0,0.5] ); cleanAxis( ax );
set( ax, 'TickDir', 'out' ); ytickangle( ax, 90 ); yticklabels( ax, yticks( ax) )
ylabel( ax, 'Density'); xlabel( ax, 'R²' )
%% R², MI, & H
mi_all = cat( 1, uMI{:} ); mi_all(isnan(mi_all))= 0;
h_all = cat(1, h{:} );
flag_groups = cat( 3,...
    [r_max_all ~= 0, r_max_all == 0], ...
 [mi_all > 0, mi_all <= 0], ...
 [h_all(:,1), ~h_all(:,1)] );
grp_tot = zeros( 8, 1 );
grp_prc = zeros( 8, 3 );
ci = 1;
logLabels = {'R', 'M', 'H'}; signflags = ones(1,3);
finalLabels = cell( 8 ,1 );
for a = flag_groups(:,1:2,1)
    for b = flag_groups(:,1:2,2)
        for c = flag_groups(:,1:2,3)
            grp_tot(ci) = sum( a & b & c );
            grp_prc(ci,:) = grp_tot(ci) ./ sum([a,b,c],1);
            finalLabels{ci} = join(logLabels);
            ci = ci + 1;
            logLabels{3} = sprintf('%c', logLabels{3} + 32*signflags(3));
            signflags(3) = -1*signflags(3);
        end
        logLabels{2} = sprintf('%c', logLabels{2} + 32*signflags(2));
        signflags(2) = -1*signflags(2);
    end
    logLabels{1} = sprintf('%c', logLabels{1} + 32*signflags(1));
    signflags(1) = -1*signflags(1);
end
%% Easy population comparison
r_squared_ppop2 = zeros( Nexp, 2 );
for cexp = 1:Nexp
    z_flag = r_max{cexp} == 0;
    aux_mu = mean( PSTHall{cexp}(:, time_flags, : ), [2,3] );
    aux_mdl = fitlm( zscore( aux_mu )', zscore( ai_pt{cexp} )', 'poly1' );
    r_squared_ppop2(cexp,1) = aux_mdl.Rsquared.Ordinary;
    aux_mu = mean( PSTHall{cexp}(:, time_flags, ~z_flag ), [2,3] );
    aux_mdl = fitlm( zscore( aux_mu )', zscore( ai_pt{cexp} )', 'poly1' );
    r_squared_ppop2(cexp,2) = aux_mdl.Rsquared.Ordinary;
end
%% Organising PSTHs by maximum R² and it's latency
rsMdl = fit_poly( [1,Nrs], [time_init, time_stop - slid_win_length], 1 );
rsq_tx = ( (1:Nrs)' .^ [1,0] ) * rsMdl;
[~, time_max] = max( r_squared_cat, [], 2 );
[H, Xedge, Yedge] = histcounts2( r_max_all, k*rsq_tx( time_max ), ...
    "BinWidth", [0.025, 10], "Normalization", "probability" );
lH = log1p( H ); ctop = round( max( lH(2:end,2:end), [], "all" ), 2 );
f = figure("Color", "w"); t = createtiles( f, 6, 6 ); ax(1) = nexttile( t, [5,1] );
histogram( ax(1), k*rsq_tx( time_max), Yedge, ...
    'Normalization', 'probability', 'EdgeColor', 'none', 'Orientation', 'horizontal')
ax(2) = nexttile( t, [5,5] );
histogram2( ax(2), 'XBinEdges', Xedge, 'YBinEdges', Yedge, ...
    'BinCounts', lH, 'EdgeColor', 'none', 'FaceColor', 'flat' );
colormap( ax(2), -gray + 1 ); clim( ax(2), [0,ctop] )
axis( ax(2), 'xy' )
aux_ax = nexttile( t, [1,1] ); set( aux_ax, 'Visible', 'off' )
ax(3) = nexttile( t, [1,5] );
histogram( ax(3), r_max_all, Xedge, ...
    'Normalization', 'probability', 'EdgeColor', 'none')
linkaxes( ax([1,2]), 'y' ); linkaxes( ax([2,3]), 'x' )
cleanAxis( ax );

title( t, 'Maximum R² and its latency per unit' )
grid( ax(2), 'off' ); ax(2).View = [0, 90];
axis( ax(2), [min(Xedge), max(Xedge), min(Yedge), max(Yedge)] )
set( get( ax(2), 'XAxis' ), 'Visible', 'off' )
set( get( ax(2), 'yAxis' ), 'Visible', 'off' )
set( get( ax(1), 'Xaxis' ), 'Scale', 'log' )
set( get( ax(3), 'Yaxis' ), 'Scale', 'log' )
xlabel( ax(3) , 'Max R²' )
ylabel( ax(1), 'Max R² latency [ms]' )
set( get( ax(1), 'XAxis' ), 'Visible', 'off' )
set( get( ax(3), 'YAxis' ), 'Visible', 'off' )
%%
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
