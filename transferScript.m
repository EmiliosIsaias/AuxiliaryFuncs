%% Roller speed per trial
rwN = 4;
[trialSubs, chTrSbs] = equiliseTrials(delayFlags, excFlag);
Mxrs = max(cellfun(@(x, y) max(x(y)), mvpt, chTrSbs'));
clrMap = bwg_cm(128);
cSb = find(contains(consCondNames,'Control'));
clMapAux = lines(Nccond);
if cSb ~= 1
    clMap = zeros(Nccond, 3); clMap(cSb,:) = 0.25*ones(1,3);
    clMap(setdiff(1:Nccond,cSb),:) = clMapAux(2:end,:);
else
    clMap = clMapAux;
end
behFig = figure('Color', 'w', 'Name', 'Behaviour');
axs = gobjects(Nccond+1,1);
axSbs = (0:rwN-2)';
for cc = 1:Nccond
    axs(cc) = subplot(rwN, Nccond, (axSbs.^[1,0])*[Nccond;cc], ...
        'Parent', behFig);
    imagesc(axs(cc), behTx*1e3, [], squeeze(vStack(:,:,trialSubs(:,cc))/Mxrs)')
    xlabel(axs(cc), 'Time [ms]')
end
ylabel(axs(1), 'Trials')
axs(Nccond+1) = subplot(rwN, Nccond, 1+Nccond*(rwN - 1):Nccond*rwN);
arrayfun(@(x) patch(axs(end), 1e3*behTx([1:end, end:-1:1]),...
    mat2ptch(rsSgnls{x}), 1, phOpts{:}, clMap(x,:)), 1:Nccond); hold on
lObj = arrayfun(@(x) plot(axs(end), 1e3*behTx, rsSgnls{x}(:,1), ...
    "Color", clMap(x,:), "LineWidth", 1.5, "DisplayName", ...
    consCondNames{x}), 1:Nccond);
xlabel(axs(end), "Time [ms]"); xlim(axs(end), 1e3*bvWin);
ylabel(axs(end), "Roller speed [cm/s]")
set(axs, axOpts{:}); title(axs(end), "Roller speed for all conditions")
lgnd = legend(axs(end), lObj); set(lgnd, lgOpts{:})
arrayfun(@(x) set(x, 'CLim', [-1,1], 'Colormap', clrMap), ...
    axs(1:end-1))
arrayfun(@(x, y) title(x, y{1}), axs(1:end-1), consCondNames')
arrayfun(@(x) set(get(x, 'YAxis'), 'Visible','off'), ...
    axs(setdiff(1:(Nccond + 1), [1, Nccond + 1])));
cb = colorbar(axs(Nccond), axOpts{1:2}, 'Location', 'west', ...
    'TickDirection', 'none', 'Ticks', [-0.9,0.9], 'TickLabels', ...
    {'Backward', 'Forward'}, 'AxisLocation', 'out');
rsPttrn = "Roller speed per trial %sVW%.2f - %.2f ms EX%s%s";
rsptName = sprintf(rsPttrn, sprintf('%s ', consCondNames{:}), bvWin*1e3, ...
    sprintf('%d ', Nex), thrshStr); rsptFile = fullfile(figureDir, rsptName);
saveFigure(behFig, rsptFile, 1)
%%
axSbs = (0:rwN-2)';
clrMap = rocket(128);
fPSTH = PSTH ./ reshape(Na*binSz, 1, 1, Nccond);
zPSTH = zscore(fPSTH, 1, [2,3]); zPopPSTH = squeeze(mean(zPSTH, 1));
Mxe = max(zPopPSTH, [], "all"); Mne = min(zPopPSTH, [], "all");
ephysFig = figure('Color', 'w', 'Name', 'Ephys');
axs = gobjects(Nccond+1,1);
for cc = 1:Nccond
    axs(cc) = subplot(rwN, Nccond, (axSbs.^[1,0])*[Nccond;cc], ...
        'Parent', ephysFig);
    imagesc(axs(cc), timeLapse*1e3, [], zPSTH(:,:,cc), [Mne, Mxe])
    xlabel(axs(cc), 'Time [ms]'); yticks(axs(cc), []);
    title(axs(cc), consCondNames{cc}); colormap(rocket(128))
end
ylabel(axs(1), 'Units');
axs(Nccond+1) = subplot(rwN, Nccond, 1+Nccond*(rwN - 1):Nccond*rwN, ...
    'NextPlot', 'add');
lObj = arrayfun(@(x) plot(axs(end), 1e3*psthTx, zPopPSTH(:,x), ...
    "Color", clMap(x,:), "LineWidth", 1.5, ...
    "DisplayName", consCondNames{x}), 1:Nccond);
lgnd = legend(axs(end), lObj); set(lgnd, lgOpts{:})
xlabel(axs(end), 'Time [ms]'); ylabel(axs(end), "Z-score")
arrayfun(@(x) set(get(x, 'YAxis'), 'Visible','off'), ...
    axs(setdiff(1:(Nccond + 1), [1, Nccond + 1])));
set(axs, axOpts{:});
cb = colorbar(axs(Nccond), axOpts{1:2}, 'Location', 'west', ...
    'TickDirection', 'none', 'AxisLocation', 'in', ...
    'Color', 0.85*ones(1,3)); cb.Label.String = 'Z-score';
ephysPttrn = 'Z-score all-units PSTH %s VW%.2f - %.2f ms Ntrials%s';
ephysName = sprintf(ephysPttrn, sprintf('%s ', consCondNames{:}), ...
    timeLapse*1e3, sprintf(' %d', Na));
ephysFile = fullfile(figureDir, ephysName);
saveFigure(ephysFig, ephysFile, 1); clearvars ephys* axs
%% Log inset for ephys and behaviour
combFig = figure('Color', 'w', 'Name','Behaviour + Ephys');
axs(1) = subplot(2, 1, 1, "Parent", combFig, "NextPlot", "add");
tmBinWidth = diff(10.^logPSTH.Log10BinEdges(:));
flPSTH = squeeze(mean(logPSTH.LogPSTH, 1))./tmBinWidth;
lObj = arrayfun(@(x) semilogx(axs(1), 1e3*logPSTH.TimeAxis, flPSTH(:,x), ...
    "Color", clMap(x,:), "LineWidth", 1.5, ...
    "DisplayName", consCondNames{x}), 1:Nccond);
lgnd = legend(lObj); set(lgnd, lgOpts{:});
title("SC Population activity in logarithmic time scale")
ylabel("Firing rate [Hz]")

axs(2) = subplot(2, 1, 2, "Parent", combFig, "NextPlot", "add");
strt = 50*1e-3;
lObj = arrayfun(@(x) semilogx(axs(2), 1e3*behTx(behTx>strt), ...
    rsSgnls{x}(behTx>strt,1), "Color", clMap(x,:), "LineWidth", 1.5, ...
    "DisplayName", consCondNames{x}), 1:Nccond);
lgnd = legend(lObj); set(lgnd, lgOpts{:});
xlabel(axs(2), "Log_{10} Time [ms]"); ylabel(axs(2), "Roller speed [cm/s]")
set(axs, axOpts{:});
arrayfun(@(x) set(get(x, 'XAxis'), 'Scale', 'log'), axs)
combPttrn = "PSTH + Roller speed%s ephysVW%.2f - %.2f behVW%.2f - %.2f";
combName = sprintf(combPttrn, sprintf(" %s", consCondNames{:}), ...
    responseWindow*1e3, strt*1e3, bvWin(2)*1e3);
combFile = fullfile(figureDir, combName);
saveFigure(combFig, combFile, 1)
%%
figure; clrmap = lines(5);
gObj = gobjects(5, 2);
for cm = 1:size(xmice.DataTable, 3)
    x = (random( normDist, size( xmice.DataTable, [1, 2]) ) + ...
        repmat( 1:3, size( xmice.DataTable, 1 ), 1 ) )';
    y = squeeze( xmice.DataTable(:,[2,1,3],cm) )';
    gObj(cm,:) = line( x, y, ...
        'Marker', 'o', 'MarkerEdgeColor', 'none', ...
        'MarkerFaceColor', clrmap(cm,:), 'Color', clrmap(cm,:), ...
        'DisplayName', xmice.MiceNames(cm));
    hold on
end
x = repmat(1:3, 10, 1);
y = BImat(:,[2,1,3]);
boxchart(x(:),  y(:), "Notch", "on", "BoxFaceColor", 0.15*ones(1,3) )
lgObj = legend(gObj(:,1), "Location", "best", "Box", "off", ...
    "Color", "none", "AutoUpdate", "off");
ylim([0,1]); title("Batch 18"); ylabel("Behaviour index"); xticks(1:3);
xticklabels({'Continuous', 'Control', 'Frequency'})
line([1,2], 1.15*[1,1]*max(BImat(:,1:2), [], "all"), 'Color', 'k', ...
    'Marker', '|')
text(1.5, 1.15*max(BImat(:,1:2), [], "all"), ...
    sprintf( "p=%.3f", signrank( behMI3(:,1) ) ), ...
    "HorizontalAlignment", "center", "VerticalAlignment", "bottom")

line(2:3, 1.15*[1,1]*max(BImat(:,[1,3]), [], "all"), 'Color', 'k', ...
    'Marker', '|')
text(2.5, 1.15*max(BImat(:,[1,3]), [], "all"), ...
    sprintf( "p=%.3f", signrank( behMI3(:,2) ) ), ...
    "HorizontalAlignment", "center", "VerticalAlignment", "bottom")

set(gca, "Box", "off", "Color", "none")

%% LFP filtering
fs = 1e3;
fo = 50;
q = 40;
bw = (fo/(fs/2))/q;
[b,a] = iircomb(fs/fo,bw,'notch');
lfp_filtered = filtfilt(b,a,lfp);
lfp_filtered = brainwaves(lfp_filtered', 1e3, {'alpha', 0.1, 100});
lfp_filtered = lfp_filtered';
figure; line(tx1k, [lfp, lfp_filtered])
lfp_z = zscore(lfp_filtered, 0);

%% MC-iRNs
Nm = size( cherryFlag, 2);
% Black, blue, baby blue
clrMap = [zeros(1,3); 0, 102, 255; 128, 179, 255]/255;

figure; hold on; arrayfun(@(c) boxchart( c + zeros( Nm, 1), ...
    mean( x_cherry{c}, 1, "omitmissing" ), "Notch", "on", ...
    "BoxFaceColor", clrMap(c,:) ), 1:3)

arrayfun(@(x) text( x, ...
    min( yv(:, x) ) * 0.85, sprintf( "n=%d", sum( ~isnan( ...
    mean( x_cherry{x}, 1, "omitmissing") ) ) ), ...
    "HorizontalAlignment", "center", "VerticalAlignment", "cap"), 1:3 )

ylim([0,1]);
jitDist = makedist("Normal", "mu", 0, "sigma", 0.1);

xv = zeros( Nm, 3); yv = xv;
for c = 1:3
    xv(:,c) = random( jitDist, Nm, 1 ) + zeros( Nm, 1) + c;
    yv(:,c) = mean( x_cherry{c}, 1, "omitmissing" );
    scatter( xv(:,c), yv(:,c), "o", "MarkerEdgeColor", "none", ...
        "MarkerFaceColor", clrMap(c, :) );
end
line( xv', yv', 'Color', 0.75 * ones(1,3) , 'LineWidth', 1/5)

p = arrayfun(@(c) signrank( mean( ...
    getMI( x_cherry{1}, x_cherry{c} ), 1, "omitmissing" ) ), 2:3 );

arrayfun(@(x) line( [1, x], ...
    max( yv(:,[1, x]), [], "all") * [1,1] * 1.15, 'Color', 'k', ...
    'Marker', '|' ), 2:3 )

posArg = {'bottom','top'};
arrayfun(@(x) text( mean( [1, x] ), ...
    max( yv(:, [1, x] ), [], "all" ) * 1.15, sprintf( "p=%.3f", p(x-1) ) , ...
    "HorizontalAlignment", "center", "VerticalAlignment", posArg{x-1} ), ...
    2:3 )

title( "MC\rightarrowiRNs behaviour index" )
ylabel("Behaviour index"); xticks(1:3); xticklabels(consCondNames)
set( gca, "Box", "off", "Color", "none", "Clipping", "off" )

%% Stan analysis
Nc = size(x, 2);
[Ns, Nm] = size(x{1}); N = Ns*Nc*Nm;
bi = zeros(N, 1, "single");
mouse_id = uint8(bi); block_id = mouse_id; tid = mouse_id;
ci = 1;
for cc = 1:size(x,2)
    for cm = 1:size(x{cc},2)
        cism = (ci:Ns+ci-1)';
        bi(cism) = x{cc}(:,cm);
        mouse_id(cism) = cm;
        block_id(cism) = (1:Ns)';
        tid(cism) = cc;
        ci = ci + Ns;
    end
end

nanflag = isnan(bi);
%% iRNs from Bayes

figure_path = fullfile( ...
    "C:\Users\neuro\seadrive_root\Emilio U", ...
    "Shared with groups\GDrive GrohLab\Projects\00 SC", ...
    "SC Behaviour\Figures\Figure 3\Matlab figures" );

% figure_path = fullfile( ...
%     "C:\Users\jefe_\seadrive_root\Emilio U", ...
%     "Für meine Gruppen\GDrive GrohLab\Projects\00 SC", ...
%     "SC Behaviour\Figures\Figure 2\Matlab figures" );

data_path = fullfile(figure_path, "Data" );

load( fullfile( data_path, "Bayes_iRN.mat") )
load( fullfile( data_path, "iRNs_4_stan.mat") )

vw = [0, 1]; binSize = 1e-2;
histOpts = {'BinLimits', vw, 'BinWidth', binSize, ...
    'Normalization', 'probability'};
fnOpts = {'UniformOutput', false};
laserSubs = 1:Nc;

bH = arrayfun(@(c) histcounts( ( params.g(:,c) * bi_sd) + bi_mu, ...
    histOpts{:}), laserSubs, fnOpts{:});
bH = cat(1, bH{:})';

mdl_bi = fit_poly( [1, size( bH, 1)], vw + [1, -1] * (binSize/2), 1);
bi_ax = ( (1:size(bH, 1))'.^[1,0] ) * mdl_bi;

X = ones( size(bH, 1), 1) * (1:Nc);
Y = bi_ax * ones( 1, Nc);
figure; patch( Y, X, bH, bH, "FaceAlpha", 1/3, "EdgeColor", "interp" )
colormap( -roma + 1 ); view([0, 69])
set( gca, "XGrid", "on", "XMinorGrid", "on")

xlim(vw); yticks(1:Nc); yticklabels({'Control', 'Laser ON', 'Laser 40 Hz'})
set(gca, "Box", "off", "Color", "none");
xlabel("Behaviour index"); title("MC\rightarrowiRNs effect on BI")

hold on; line( median( (params.g(:,laserSubs) * bi_sd) + bi_mu ), 1:Nc,...
    'LineWidth', 2, 'Color', 'k', 'Marker', 'x', 'LineStyle', 'none' )

cb = colorbar("Box", "off", "Location", "westoutside", "Ticks", [] );
cb.Label.String = "Low \leftarrow BI likelihood \rightarrow High";

set( get( gca, "ZAxis"), "visible", "off")

%%
for a = 1:size( aH, 2 )
    figure(a);
    contour( ...
        ones(size(aH{a}, 2), 1) * ( 1:size(aH{a}, 1) ), ...
        eff_ax(:) * ones(1,size(aH{a},1)), aH{a}', 10, 'DisplayName', ...
        'Likelihood');
    colormap(-gray(10)+1); ylabel( 'Movement bias' )
    xticks( 1:size( aH{a}, 1 ) ); xticklabels( consCondNames )
    title( sprintf( "Mouse %d (%s)", a, ...
        strrep(xmice_iRNs.MiceNames(a), '_', '\_' ) ) )
    grid on; set(gca, "Box", "off", "Color", "none", "Clipping", "off")
    cb = colorbar("peer", gca, ...
        "eastoutside", "AxisLocation", "out", "Box", "off");
    hold on; line( 1:size( aH{a}, 1 ), ...
        squeeze( mean( params.alpha(:, a, :) ) ), 'Color', 'k', ...
        'LineWidth', 2 , 'DisplayName', 'Mean')
    line( 1:size( aH{a}, 1 ), ...
        squeeze( median( params.alpha(:, a, :) ) ), 'Color', 'b', ...
        'LineWidth', 1 , 'DisplayName', 'Median', 'Marker', 'x')
    lgObj = legend("show", "Box", "off", "Color", "none");
end

%% All mice
figurePath = ...
    fullfile( "C:\Users\jefe_\seadrive_root\Emilio U\Für meine Gruppen\GDrive GrohLab\Projects\00 SC\SC Behaviour\Figures\Figure 3\Matlab figures" );

c_MI = arrayfun(@(c) getMI( params.alpha(:,:,1), params.alpha(:, :, c) ), ...
    2:3, fnOpts{:} ); c_MI = cat(3, c_MI{:} );

vw = [-1,1];
histOpts = {'BinLimits', vw, 'BinWidth', binSize, 'Normalization', 'probability'};
mH = arrayfun(@(a) arrayfun(@(c) histcounts( c_MI(:, a, c), ...
    histOpts{:}), 1:2, fnOpts{:}), 1:size(aH, 2), fnOpts{:});
mH = cellfun(@(cc) cat( 1, cc{:} ), mH, fnOpts{:});
mH = cat( 3, mH{:} );

mdlMI = fit_poly( [1, size(mH,2)], vw + [1,-1]*(binSize/2), 1);
mi_ax = ( (1:size(mH, 2))'.^[1, 0] ) * mdlMI;
X =  ones( size(mH, 2), 1 ) * ( 1:size(mH,3) );
Y = mi_ax(:) * ones(1, size(mH, 3) );
figs = gobjects(2, 1);
for cmi = 1:2
    figs(cmi) = figure('Name', sprintf( "MI for %s", consCondNames{cmi+1} ) );
    contour( X, Y, squeeze( mH(cmi, :, :) ), 10);
    colormap( viridis );
    yline( 0 , 'r'); ylim( vw ); ylabel( 'Modulation index' )
    xticks( 1:size( aH, 2) )
    xticklabels( strrep( xmice_iRNs.MiceNames, '_', '\_' ) )
    title( sprintf( 'Modulation index for %s w.r.t. control', ...
        consCondNames{cmi+1}) )
    cb = colorbar("Box", "off", "Location", "eastoutside");
    cb.Label.String = "MI likelihood";
    set( gca, "Box", "off", "Color", "none", "Clipping", "off" )
    hold on; line( 1:size( aH, 2 ), squeeze( median( c_MI( :, :, cmi ), 1 ) ), ...
        'Color', 'k', 'LineStyle', 'none' , 'DisplayName', 'Median', ...
        'Marker', 'x', 'LineWidth', 2);
    saveFigure(figs(cmi), fullfile( figurePath, ...
        sprintf("MI_%s_Bayes", consCondNames{cmi+1} ) ), true )
end

%% Lick bayes

load( "C:\Users\neuro\seadrive_root\Emilio U\Shared with groups\GDrive GrohLab\Projects\00 Salience\Bayes Model Data\lick_tables_Native_DREADDs_IR_RS_all_sessions_with_trialID.mat" )
load( "Z:\Nadine\Behavior_Analysis\lick_tables_Native_DREADDs_IR_RS_all_sessions_with_trialID.mat" )

fnOpts = {'UniformOutput', false};
getAllTrialsFromMice = @(tbl, m) tbl{ string( tbl{:, "MouseID"} ) == ...
    string( m ), "TrialID" };
getMaxMiceTrialNum = @(tbl, m) max( getAllTrialsFromMice(tbl, m) );
lickOrNotToLick = @(tbl) isnan( tbl{:, "Lick_Latencies"});
decodeTreatment = @(tbl) uint8( contains( tbl{:, "Label"}, {'Miss', 'Hit'}) + 1 );
getLickTime = @(tbl) tbl{:, "Lick_Latencies"};

nmice_original_id = unique( Native_IR_lick_table{:,"MouseID"} );
[~, nir_mouse_id] = ismember( Native_IR_lick_table{:, "MouseID"}, ...
    nmice_original_id ); nir_mouse_id = uint8( nir_mouse_id );
[~, nrs_mouse_id] = ismember( Native_RS_lick_table{:, "MouseID"}, ...
    nmice_original_id ); nrs_mouse_id = uint8( nrs_mouse_id );

Nt_nir = uint16( cellfun(@(m) getMaxMiceTrialNum( Native_IR_lick_table, m), ...
    nmice_original_id ) );
Nt_nir_progress = arrayfun(@(m) ...
    getAllTrialsFromMice( Native_IR_lick_table, nmice_original_id{m} )./ ...
    single( Nt_nir(m) ), ( 1:size(nmice_original_id, 1) )' , fnOpts{:});
Nt_nir_progress = cat(1, Nt_nir_progress{:} );

Nt_nrs = uint16( cellfun(@(m) getMaxMiceTrialNum( Native_RS_lick_table, m), ...
    nmice_original_id ) );
Nt_nrs_progress = arrayfun(@(m) ...
    getAllTrialsFromMice( Native_RS_lick_table, nmice_original_id{m} )./ ...
    single( Nt_nrs(m) ), ( 1:size(nmice_original_id, 1) )' , fnOpts{:});
Nt_nrs_progress = cat(1, Nt_nrs_progress{:} );


dmice_original_id = unique( DREADDs_IR_lick_table{:, "MouseID"} );
[~, dir_mouse_id] = ismember( DREADDs_IR_lick_table{:, "MouseID"}, ...
    dmice_original_id ); dir_mouse_id = uint8( dir_mouse_id ) + ...
    size( nmice_original_id, 1);
[~, drs_mouse_id] = ismember( DREADDs_RS_lick_table{:, "MouseID"}, ...
    dmice_original_id ); drs_mouse_id = uint8( drs_mouse_id ) + ...
    size( nmice_original_id, 1);

Nt_dir = uint16( cellfun(@(m) getMaxMiceTrialNum( DREADDs_IR_lick_table, m ), ...
    dmice_original_id ) );
Nt_dir_progress = arrayfun(@(m) ...
    getAllTrialsFromMice( DREADDs_IR_lick_table, dmice_original_id{m} )./ ...
    single( Nt_dir(m) ), ( 1:size(dmice_original_id, 1) )' , fnOpts{:});
Nt_dir_progress = cat(1, Nt_dir_progress{:} );

Nt_drs = uint16( cellfun(@(m) getMaxMiceTrialNum( DREADDs_RS_lick_table, m ), ...
    dmice_original_id ) );
Nt_drs_progress = arrayfun(@(m) ...
    getAllTrialsFromMice( DREADDs_RS_lick_table, dmice_original_id{m} )./ ...
    single( Nt_drs(m) ), ( 1:size(dmice_original_id, 1) )' , fnOpts{:});
Nt_drs_progress = cat(1, Nt_drs_progress{:} );


mouse_id = cat(1, nir_mouse_id, nrs_mouse_id, dir_mouse_id, drs_mouse_id );
progress_id = uint16( round (cat(1, Nt_nir_progress*1e3, ...
    Nt_nrs_progress*1e3, Nt_dir_progress*1e3, Nt_drs_progress*1e3 ) ) );
lick_flag = cat(1, lickOrNotToLick(Native_IR_lick_table), ...
    lickOrNotToLick(Native_RS_lick_table), ...
    lickOrNotToLick(DREADDs_IR_lick_table), ...
    lickOrNotToLick(DREADDs_RS_lick_table) );
tid = cat(1, decodeTreatment(Native_IR_lick_table), ...
    decodeTreatment(Native_RS_lick_table), ...
    decodeTreatment(DREADDs_IR_lick_table), ...
    decodeTreatment(DREADDs_RS_lick_table) );
block_id = [ones( sum( Nt_nir ), 1, "uint8"); ...
    ones( sum( Nt_nrs ), 1, "uint8")+1; ...
    ones( sum( Nt_dir ), 1, "uint8"); ...
    ones( sum( Nt_drs ), 1, "uint8")+1];

lick_lat = cat(1, getLickTime(Native_IR_lick_table), ...
    getLickTime(Native_RS_lick_table), ...
    getLickTime(DREADDs_IR_lick_table), ...
    getLickTime(DREADDs_RS_lick_table) );

N = uint16( length( lick_flag ) ); Nts = uint8( max( tid ) );
Nm = uint8( max( mouse_id ) );
%{
save( fullfile( "C:\Users\jefe_\seadrive_root\Emilio U\Für meine Gruppen", ...
    "GDrive GrohLab\Projects\00 Salience\Bayes Model Data", ...
    "Lick_data_4_R.mat"), "mouse_id", "block_id", "tid", "lick_flag", ...
    "progress_id", "N", "Nm", "Nts", "lick_lat" ) 
%}
save( "Z:\Nadine\Behavior_Analysis\Bayes\Lick_Bayes.mat", ...
    "mouse_id", "block_id", "tid", "lick_flag", ...
    "progress_id", "N", "Nm", "Nts", "lick_lat" )

%% Pharmacology effects from Bayes
%{
vw = [-2, 2]; binSize = 1e-2;
histOpts = {'BinLimits', vw, 'BinWidth', binSize, ...
    'Normalization', 'probability'};
fnOpts = {'UniformOutput', false};

gH = arrayfun(@(c) histcounts( params.g(:,c), histOpts{:}), 1:3, fnOpts{:});
gH = cat(1, gH{:})';

mdl_zbi = fit_poly( [1, size( gH, 1)], vw + [1, -1] * (binSize/2), 1);
zbi_ax = ( (1:size(gH, 1) )'.^[1,0] ) * mdl_zbi;

X = ones( size(gH, 1), 1) * (1:3);
Y = zbi_ax * ones( 1, 3);
figure; contour( X, Y, gH, 10); colormap( viridis )
%}

figure_path = fullfile( ...
    "C:\Users\neuro\seadrive_root\Emilio U", ...
    "Shared with groups\GDrive GrohLab\Projects\00 SC", ...
    "SC Behaviour\Figures\Figure 1\Matlab figures" );
data_path = fullfile( figure_path, "Data");

load( fullfile( data_path, "Pharma.mat" ))

vw = [0, 1]; binSize = 1e-2;
histOpts = {'BinLimits', vw, 'BinWidth', binSize, ...
    'Normalization', 'probability'};
fnOpts = {'UniformOutput', false};
pharmaSubs = [2,1,3];
consCondNames = {'Muscimol', 'Control', 'Picrotoxin'};
bH = arrayfun(@(c) histcounts( ( params.g(:,c) * bi_scale) + bi_centre, ...
    histOpts{:}), pharmaSubs, fnOpts{:});
bH = cat(1, bH{:})';

mdl_bi = fit_poly( [1, size( bH, 1)], vw + [1, -1] * (binSize/2), 1);
bi_ax = ( (1:size(bH, 1))'.^[1,0] ) * mdl_bi;

X = ones( size(bH, 1), 1) * (1:3);
Y = bi_ax * ones( 1, 3);
figure; patch( Y, X, bH, bH, "FaceAlpha", 1/3, "EdgeColor", "interp" )
xlim(vw); yticks(1:3); yticklabels({'Muscimol', 'Control', 'Picrotoxin'})
colormap( -roma + 1); view([0, 60]);
set( gca, "XGrid", "on", "XMinorGrid", "on")

set(gca, "Box", "off", "Color", "none");
xlabel("Behaviour index"); title("Effects of silencing / desinhibiting SC")

hold on; line( median( (params.g(:,pharmaSubs) * bi_scale) + bi_centre ), 1:3,...
    'LineWidth', 2, 'Color', 'k', 'Marker', 'x', 'LineStyle', 'none' )

cb = colorbar("Box", "off", "Location", "westoutside", "Ticks", [] );
cb.Label.String = "Low \leftarrow BI likelihood \rightarrow High";
set( get( gca, "ZAxis"), "visible", "off")

%% eOPN3 & ChR2
figure_path = fullfile( ...
    "C:\Users\neuro\seadrive_root\Emilio U", ...
    "Shared with groups\GDrive GrohLab\Projects\00 SC", ...
    "SC Behaviour\Figures\Figure 2\Matlab figures" );

% figure_path = fullfile( ...
%     "C:\Users\jefe_\seadrive_root\Emilio U", ...
%     "Für meine Gruppen\GDrive GrohLab\Projects\00 SC", ...
%     "SC Behaviour\Figures\Figure 2\Matlab figures" );

data_path = fullfile(figure_path, "Data" );

load( fullfile( data_path, 'eOPN3_pool.mat' ) )
load( fullfile( data_path, 'ChR2_pool.mat' ) )

Nmo = size( oPN3, 1); Nmc = size( chr2, 1);
Nc = size( oPN3, 2 ) + size( chr2, 2 ) - 2;

om_ids = 1:6; cm_ids = 7:18;

bi = [reshape( oPN3(:, 1:2), [], 1); reshape( chr2, [], 1)];
[~, bi_centre, bi_scale] = zscore(bi(~isnan(bi)), 0);

omouse_id = reshape( ( om_ids )' * ones(1, 2), [], 1 );
cmouse_id = reshape( ( cm_ids )' * ones(1, 7), [], 1 );
mouse_id = uint8( [omouse_id; cmouse_id] );

otid = reshape( ones( Nmo, 1) * (1:2), [], 1 );
ctid = reshape( ones( Nmc, 1) * [1,3:Nc], [], 1 );
tid = uint8( [otid; ctid ]);

save( fullfile( data_path, "MC_opsines_4_r.mat" ), "Nc", "mouse_id", ...
    "tid", "bi", "bi_centre", "bi_scale")

%% MC-Opsines Bayes
figure_path = fullfile( ...
    "C:\Users\neuro\seadrive_root\Emilio U", ...
    "Shared with groups\GDrive GrohLab\Projects\00 SC", ...
    "SC Behaviour\Figures\Figure 2\Matlab figures" );

% figure_path = fullfile( ...
%     "C:\Users\jefe_\seadrive_root\Emilio U", ...
%     "Für meine Gruppen\GDrive GrohLab\Projects\00 SC", ...
%     "SC Behaviour\Figures\Figure 2\Matlab figures" );

data_path = fullfile(figure_path, "Data" );

load( fullfile( data_path, "Opsines_Bayes.mat") )
%%
vw = [0, 1]; binSize = 1e-2;
histOpts = {'BinLimits', vw, 'BinWidth', binSize, ...
    'Normalization', 'probability'};
fnOpts = {'UniformOutput', false};
opsinSubs = [2,1,3:size(params.g,2)];

bH = arrayfun(@(c) histcounts( ( params.g(:,c) * bi_scale) + bi_centre, ...
    histOpts{:}), opsinSubs, fnOpts{:});
bH = cat(1, bH{:})';

mdl_bi = fit_poly( [1, size( bH, 1)], vw + [1, -1] * (binSize/2), 1);
bi_ax = ( (1:size(bH, 1))'.^[1,0] ) * mdl_bi;

X = ones( size(bH, 1), 1) * (1:Nc);
Y = bi_ax * ones( 1, Nc);
figure; patch( Y, X, bH, bH, "FaceAlpha", 1/3, "EdgeColor", "interp" )
colormap( -roma + 1); view([0, 69])
xlim(vw); yticks(1:Nc); yticklabels({'eOPN3', 'Control', '30 ms', ...
    '30 ms 40 Hz', '100 ms', '100 ms 40 Hz', '400 ms', '400 ms 40 Hz'})
set(gca, "Box", "off", "Color", "none");
xlabel("Behaviour index"); title("Effects of inhibiting / activating MC")

hold on; line( median( (params.g(:,opsinSubs) * bi_scale) + bi_centre ), 1:Nc,...
    'LineWidth', 2, 'Color', 'k', 'Marker', 'x', 'LineStyle', 'none' )

cb = colorbar("Box", "off", "Location", "westoutside", "Ticks", [] );
cb.Label.String = "Low \leftarrow BI likelihood \rightarrow High";

set( get( gca, "ZAxis"), "visible", "off")


%% Z-maximum value
Na = sum(behData.Conditions);
normDist = makedist("Normal", "mu", 0, "sigma", 1 );
z_mvpt = arrayfun(@(x) cat(1, x.Results.Z_Amplitude)', behRes, fnOpts{:});
z_mvpt = cat(1, z_mvpt{:});
ylmts = [min( z_mvpt(:) ), max( z_mvpt(:) )] * 1.1;
zax = ylmts(1):range(z_mvpt(:))/100:ylmts(2);
figure("Color", "w"); %imagesc( [1,4], ylmts, pdf( normDist, zax )', ...
%"AlphaData", 1/3);
hold on; colormap(-gray + 1); xlim([0,5])
bxObj = boxchart( reshape( ones(sum( Na ), 1) * (1:Nb), [], 1), ...
    reshape( z_mvpt, [], 1), "Notch", "on", ...
    "GroupByColor", repmat( behData.Conditions, 4, 1 ) * (1:Nccond)', ...
    "JitterOutliers", "on", "MarkerStyle", "." );
xticks(1:4); xticklabels(behNames);
set( gca, "Box", "off", "Color", "none" );
legend(bxObj, consCondNames, "Color", "none", "Box", "off", "Location", "best")
title("Z-score of maximum per trial"); ylabel("Z-Max")
%% Find thresholds for a z-distribution
alph = 0:0.01:1;
sigmaTh = arrayfun(@(z) ...
    fminbnd(@(y) ...
    norm( integral(@(x) normDist.pdf(x), -y, y) - z, 2 ), 0, 5 ), ...
    alph );
% figure; plot( alph, sigmaTh )

z_flags = z_mvpt > reshape( -[sigmaTh, inf], 1, 1, [] ) & ...
    z_mvpt < reshape( [sigmaTh, inf], 1, 1, [] );

z_th_curves = arrayfun(@(c) ...
    squeeze( sum( z_flags( pairedStimFlags(:,c), :, :) ) )./Na(c), ...
    1:Nccond, fnOpts{:} );
figu
z_auc = cellfun(@(z) sum( z, 2 )./ (numel( sigmaTh ) + 1), ...
    z_th_curves, fnOpts{:} );

figure("Color", "w");
for cbp = 1:4
    ax = subplot( 2, 2, cbp, "NextPlot", "add" );
    arrayfun(@(c) line( ax, [sigmaTh, sigmaTh(end)+1] , ...
        z_th_curves{c}(cbp,:)' ), 1:Nccond)
    title(ax, behNames(cbp));
    if cbp==1
        ylabel(ax, 'Trial proportion')
    elseif cbp==4
        xlabel(ax, '\Theta_Z')
    end
    set(ax, axOpts{:})
    legend( arrayfun(@(c) ...
        join( [consCondNames(c), string( z_auc{c}(cbp) )]), 1:Nccond ), ...
        lgOpts{:})
end
%% Pooling script (Muscimol and PTX)
% Muscimol mice
load('Z:\Emilio\SuperiorColliculusExperiments\Roller\Batch6_beh+Muscimol\Batch6_BehaviourIndex.mat')
mice6 = mice;
load('Z:\Emilio\SuperiorColliculusExperiments\Roller\Batch11_ephys.MC\Batch11_BehaviourIndex.mat')
mice11 = mice;

% PTX mice
load('Z:\Emilio\SuperiorColliculusExperiments\Roller\Batch12_ephys.e\Batch12_BehaviourIndex.mat')
mice12 = mice;
load('Z:\Emilio\SuperiorColliculusExperiments\Roller\Batch13_beh\Batch13_BehaviourIndex.mat')
mice13 = mice;
load('Z:\Emilio\SuperiorColliculusExperiments\Roller\Batch10_ephys.e\Batch10_BehaviourIndex.mat')
mice10 = mice;

dom = -5:0.05:5;
mpdf = @(x) pdf( x, dom );
getSessPerBatch = @(ms) arrayfun(@(m) numel(m.Sessions) , ms );
unpackTbl = @(tbl) table( tbl.Conditions{:}, tbl.Trial_and_Amp_Indices{:}, ...
    tbl.PolygonUnfold{:}, 'VariableNames', tbl.Properties.VariableNames );
getCondsInBatch = @(mc) arrayfun(@(m) arrayfun(@(s) ...
    size( s.DataTable, 1), m.Sessions), mc );
total_var_dist = @(dmat) integral( @(x) ...
    abs( pdf( dmat(1), x ) - pdf( dmat(2), x ) ), dom(1), dom(end) );
kl_div = @(dmat) KullbackLeiblerDivergence( mpdf( dmat(1) ), mpdf( dmat(2) ) );
fnOpts = {'UniformOutput', false}; 


mice_mus = [mice6; mice11]; Nmm = numel( mice_mus );
mice_ptx = [mice10; mice12; mice13]; Nmp = numel( mice_ptx );


Ncondm = getCondsInBatch( mice_mus );
musAreas = zeros( max(Ncondm), Nmm, 2 );
for cmm = 1:Nmm
    for cs = 1:numel(mice_mus(cmm).Sessions)
        musAreas(1:Ncondm, cmm, :) = reshape( ...
            mice_mus(cmm).Sessions(cs).DataTable.Trial_and_Amp_Indices, ...
            Ncondm(cmm), 1 , 2 );
    end
end

Ncondp = getCondsInBatch( mice_ptx );
ptxAreas = nan( max(Ncondp), Nmp, 2 );
for cmm = 1:Nmp
    for cs = 1:numel(mice_ptx(cmm).Sessions)
        aux_table = mice_ptx(cmm).Sessions(cs).DataTable;
        ptxAreas(1:Ncondp(cmm), cmm, :) = reshape( ...
            aux_table.Trial_and_Amp_Indices, Ncondp(cmm), 1 , 2 );
    end
end
%%
tvp_mus = arrayfun(@(bp) arrayfun(@(m) total_var_dist( bDist_mus(:, bp, m) ), 1:Nmm), 1:Nbs, fnOpts{:} );
tvp_mus = cat( 1, tvp_mus{:} );
kl_mus = arrayfun(@(bp) arrayfun(@(m) kl_div( bDist_mus(:, bp, m) ), 1:Nmm), 1:Nbs, fnOpts{:} );
kl_mus = cat( 1, kl_mus{:} );
%%
Nconds_mus = arrayfun(@(m) arrayfun(@(s) size( s.DataTable, 1 ), m.Sessions), mice_mus );
Nconds_ptx = arrayfun(@(m) arrayfun(@(s) size( s.DataTable, 1 ), m.Sessions), mice_ptx );
prmSubs_mus = arrayfun(@(x) nchoosek(1:x,2), Nconds_mus, fnOpts{:} );
prmSubs_ptx = arrayfun(@(x) nchoosek(1:x,2), Nconds_ptx, fnOpts{:} );
Ncombs_mus = cellfun(@(pr) size( pr, 1 ), prmSubs_mus );
Ncombs_ptx = cellfun(@(pr) size( pr, 1 ), prmSubs_ptx );

Nbs = cellfun(@(d) size( d, 2), bDist );
tvd = cell( Nbr, 1);
for cc = 1:Nbr
    if Nconds(cc) > 1
        tvd{cc} = zeros( Nbs(cc), Ncombs(cc) );
        for cr = 1:Ncombs(cc)
            ps = prmSubs{cc}(cr,:);
            for cbp = 1:Nbs(cc)
                tvd{cc}(cbp, cr) = total_var_dist( bDist{cc}(cbp, ps) );
            end
        end
    else
        fprintf(1, 'Unsure what to do\n')
    end
end

ptx_vals = arrayfun(@(m) cellfun(@(p) permute( p', 3:-1:1 ), ...
    m.Sessions.DataTable.PolygonUnfold', fnOpts{:} ), mice_ptx, fnOpts{:} );
p_p = arrayfun(@(c) signrank( ptx_vals2{1}(:,c,1), ptx_vals2{2}(:,c,1) ), 1:8 );
ptx_vals2 = cell( size( ptx_vals, 1), 3 );
for cr = 1:size( ptx_vals, 1)
    ptx_vals2(cr,1:size(ptx_vals{cr},2)) = ptx_vals{cr,:};
end
ptx_vals = ptx_vals2;
ptx_vals2 = arrayfun(@(c) cat(1, ptx_vals{:,c} ), 1:2, fnOpts{:} );
p_p = arrayfun(@(c) signrank( ptx_vals2{1}(:,c,1), ptx_vals2{2}(:,c,1) ), 1:8 );

bp_idx = tocol( repmat( ones( Nmp , 1) * (1:8), 1, 2 ) );
treat_idx = tocol( ones( 8*Nmp, 1)*(1:2) );

bxObj = boxchart(ax, bp_idx, bpa_ptx_vec, 'Notch', 'on', 'JitterOutliers', 'on', 'GroupByColor', treat_idx, 'MarkerStyle', '.' );
xticks( 1:8 )
xlim(ax, [0.5, 8.5] )
xticklabels( ax, bodypart_names )
bxObj(2).BoxFaceColor = ones(1, 3);
bxObj(2).BoxEdgeColor = zeros(1, 3);
bxObj(1).BoxFaceColor = zeros(1,3 );
%%
figs = gobjects( 8, 1 );
newFig = @(n) figure('Color', 'w', 'Name', n);
tiles = @(x,r,c) tiledlayout(x, r, c, 'TileSpacing','tight','Padding', 'tight');
lnOpts = {'LineStyle', 'none', 'Marker', 'o', ...
    'MarkerFaceColor', 0.35*ones(1,3), 'MarkerEdgeColor', 'none'};
catCols = @(x) cat(2, x{:});
fetchValues = @(m, b, s) catCols( cellfun(@(x) x(:,b,s), m, fnOpts{:} ) );
my_scatt = @(aax, mat) line( aax, mat(:,1), mat(:,2), lnOpts{:});
yeqxLine = @(x) line(x, xlim(x), xlim(x), 'LineStyle', '--', ...
            'Color', 0.15*ones( 1, 3 ) );
yLabels = {'Trial proportion', 'Amplitude index'};
drug_names = {'Muscimol', 'Picrotoxin'};
cleanAxis = @(x) set( x, "Box", "off", "Color", "none" );
for cf = 1:8
    figs(cf) = newFig( cf, bodypart_names(cf) );
    t = tiles( figs(cf) );
    for cc = 1:4
        ax = cleanAxis( nexttile(t) ); ax.NextPlot = 'add';
        cm = ceil( cc/2 );
        if cm == 1
            title( ax, drug_names{cc} )
        end
        if mod(cc,2)
            aux_mat = fetchValues( musc_vals, cf, cm );
            ylabel(ax, join( [string( yLabels{cm} ), "drug"] ) )
        else
            aux_mat = fetchValues( ptx_vals2, cf, cm );
            if cm == 2
                xlabel(ax, "Control" )
            end
        end
        my_scatt(ax, aux_mat); yeqxLine(ax)
    end
    title( t, bodypart_names(cf) )
end

%%
dotOutliers = @(x) set( x, 'MarkerStyle', '.', 'MarkerColor', 'k' );
fig = newFig( "Baseline PTX" ); t = tiles( fig, 2, 1 );
ax = nexttile( t ); cleanAxis( ax );
bp_idx = ones( size( bDist_ptx, 1), 1 ) * (1:Nbp );
bp_idx = bp_idx(:);
treat_idx = ones( numel( bp_idx ), 1) * (1:2);
treat_idx = treat_idx(:);
bxObj = boxchart( ax, repmat( bp_idx, 2, 1), disp_ptx(:), 'GroupByColor', treat_idx, 'Notch', 'on', 'JitterOutliers', 'on' );
bxObj(1).BoxFaceColor = 'k';
bxObj(2).BoxFaceColor = 'none';
bxObj(2).BoxEdgeColor = zeros(1, 3);
arrayfun(dotOutliers, bxObj );
ylabel(ax, 'Norm. IQR')
ax.XAxis.Visible = "off";
lgObj = legend({'Control', 'PTX'}, lgOpts{:} );
ax = nexttile( t ); cleanAxis( ax );
bxObj = boxchart( ax, repmat( bp_idx, 2, 1), med_ptx(:), 'GroupByColor', treat_idx, 'Notch', 'on', 'JitterOutliers', 'on' );
ylabel(ax, 'Norm. median')
bxObj(1).BoxFaceColor = 'k';
bxObj(2).BoxFaceColor = 'none';
bxObj(2).BoxEdgeColor = zeros(1, 3);
arrayfun(dotOutliers, bxObj );
xticks(ax, 1:Nbp )
xticklabels(ax, bodypart_names )
%%
roller_path = "Z:\Emilio\SuperiorColliculusExperiments\Roller";
bp_paths = dir( fullfile( roller_path, "\Batch*\Batch*_BehaviourIndex.mat" ) );
bodypart_names = string( {bp_paths.name} )';
b_num = regexp( bodypart_names, '\d+', 'match' );
b_num = cat(1, b_num{:} ); b_num = str2double( b_num );
expandPath = @(x) fullfile( x.folder, x.name);
fnOpts = {'UniformOutput', false};

mc_subs = [1, 2, 7, 8, 10, 12, 15, 16, 18];
bc_subs = [1, 2, 19];
bs_subs = 2;
mt_subs = [11, 14, 17];

exp_type_subs_cell = {mc_subs, bc_subs, bs_subs, mt_subs};
clearvars *_subs -except exp_type_subs_cell

mice_bulk = arrayfun( @(x) load( expandPath( x ), "mice" ), ...
    bp_paths );
mice_exp_sub = cellfun(@(x) any(b_num == x, 2), exp_type_subs_cell, ...
    fnOpts{:} );

summMice = cellfun(@(x) summariseMiceBeh( cat(1, mice_bulk(x).mice ) ), ...
    mice_exp_sub, fnOpts{:} );
