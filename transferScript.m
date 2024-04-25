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
figurePath = fullfile( "C:\Users\jefe_\seadrive_root\Emilio U\Für meine Gruppen\GDrive GrohLab\Projects\00 SC\SC Behaviour\Figures\Figure 3\Matlab figures" );

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

%% Analysing behaviour alone

exp_path = fullfile("Z:\Emilio\SuperiorColliculusExperiments\Roller\Batch13_beh\PTX\WT2\230324_PiTX");
expandPath = @(x) fullfile( x.folder, x.name);
m = 1e-3;

eph_path = dir( fullfile( exp_path, "ephys*" ) );
beh_path = fullfile( exp_path, "Behaviour" );
if ~isempty( eph_path )
    eph_path = expandPath( eph_path );
    figure_path = fullfile( eph_path, "Figures" );
    af_path = dir( fullfile( eph_path , "*analysis.mat" ) );
    af_path = expandPath( af_path );
elseif exist( beh_path, "dir" )
    af_path = expandPath( dir( fullfile( beh_path, "*analysis.mat") ) );
    figure_path = fullfile( beh_path, "Figures" );
else
    beh_path = exp_path;
    af_path = expandPath( dir( fullfile( beh_path, "*analysis.mat") ) );
    figure_path = fullfile( beh_path, "Figures" );
end

[~, af_name] = fileparts( af_path );
expName = extractBefore(af_name, "analysis");
load( af_path, "Conditions", "fs")

fnOpts = {'UniformOutput', false};
hstOpts = {'BinMethod', 'integers', 'BinLimits', [-0.5,4.5]};
axOpts = {'Box','off','Color','none'};
lgOpts = cat( 2, axOpts{1:2}, {'Location','best'} );

open Conditions
ldFlag = false;
try
    load( expandPath( dir( fullfile( beh_path, "RollerSpeed*.mat" ) ) ), "fr")
catch
    ldFlag = true;
end

%% Run independently
% User input!!
consCond = 2:3;
Nccond = length( consCond );
prmSubs = nchoosek(1:Nccond,2);

pairedStimFlags = arrayfun(@(c) any( ...
    Conditions(1).Triggers(:,1) == ...
    reshape( Conditions(c).Triggers(:,1), 1, []), 2), consCond, fnOpts{:} );
pairedStimFlags = cat(2, pairedStimFlags{:});

consCondNames = string( { Conditions( consCond ).name  } );

[behRes, behFig_path, behData, aInfo] = analyseBehaviour(beh_path, ...
    "ConditionsNames", cellstr(consCondNames), ...
    "PairedFlags", pairedStimFlags, ...
    "FigureDirectory", figure_path, ...
    "ResponseWindow", [25, 350] * m, ...
    "ViewingWindow", [-450, 500] * m);

if ~exist( "fr", "var") && ldFlag
    load( expandPath( dir( fullfile( beh_path, "RollerSpeed*.mat" ) ) ), "fr")
    ldFlag = false;
end

vwin = sscanf( aInfo.VieWin, "V%f - %f s")';
mdlt = fit_poly( [1, size( behData.Data, 1)], vwin + [1,-1] * (1/(2 * fr) ), 1);
txb = ( (1:size( behData.Data, 1))'.^[1,0] ) * mdlt;
behNames = string( { behRes(1).Results.BehSigName } );
 
biFigPttrn = "BehIndex%s";
biFigPttrn = sprintf(biFigPttrn, sprintf(" %s (%%.3f)", consCondNames));
[pAreas, ~, behAreaFig] = createBehaviourIndex(behRes);
behRes = arrayfun(@(bs, ba) setfield(bs,'BehIndex', ba), behRes, pAreas);
set(behAreaFig, 'UserData', behRes)

biFN = sprintf(biFigPttrn, pAreas);

trMvFlag = arrayfun(@(cr) behRes(1).Results(cr).MovStrucure.MovmentFlags, ...
    1:size(behRes(1).Results,2), fnOpts{:}); trMvFlag = cat(3, trMvFlag{:});
BIscaleMat = sum(trMvFlag,3);
BIscale = arrayfun(@(cc) BIscaleMat(pairedStimFlags(:,cc), cc), 1:Nccond, ...
    fnOpts{:});
[hg, hg_bin] = cellfun(@(c) histcounts(c, hstOpts{:}), ...
    BIscale, fnOpts{:});
hg = cat(1, hg{:}); hg_bin = cat(1, hg_bin{:});

% [p_amp, h_amp] = ranksum(cat(1, zamp{1,:}), cat(1, zamp{2,:}));

clrMap = lines(Nccond);
countFig = figure; ax(1) = subplot(10,1,1:8);
bar(ax(1), (0:4)', (hg./sum(hg,2))', 'EdgeColor', 'none'); hold on;
poaDist = cellfun(@(bi) fitdist(bi,"Poisson"), BIscale);
ylim(ax(1), [0,1]); set(ax(1), axOpts{:});
legend(ax(1), consCondNames, 'AutoUpdate','off', lgOpts{:})
lmbdaHeight = 0.95-(0.15/Nccond)*(0:Nccond-1);
arrayfun(@(pd) scatter(ax(1), poaDist(pd).lambda, lmbdaHeight(pd), '|',...
    'MarkerEdgeColor', clrMap(pd,:)), 1:Nccond)
arrayfun(@(pd) line(ax(1), paramci(poaDist(pd)), ...
    lmbdaHeight([pd,pd]), 'Color', clrMap(pd,:), ...
    'Marker', '|'), 1:Nccond)
[p, chiVal] = arrayfun(@(ps) chi2test(hg(prmSubs(ps,:), :)), ...
    1:size(prmSubs,1));
ax(2) = subplot(10,1,9:10);
signBeh = arrayfun(@(x) sprintf("%s vs %s p=%.3f", ...
    consCondNames(prmSubs(x,:)), p(x)), 1:size(prmSubs,1));
text(ax(2), 0, -0.3, sprintf('%s vs. %s P=%.3f\n', ...
    [consCondNames(prmSubs), string(p(:))]'))
set(ax(2), 'Visible', 'off')
set(countFig, 'UserData', {signBeh, p})
title(ax(1), strrep(expName, '_',' ')); xlabel(ax(1),'Moving body parts')
ylabel(ax(1),'Trial proportion')
countFigName = sprintf("Count distributions P%s", ...
    sprintf(" %.3f", p(:)));

saveFigure(behAreaFig, fullfile(behFig_path, biFN), true);
saveFigure(countFig, fullfile(behFig_path, countFigName), true);

%% Normalised amplitud by absolute maximum

fig = figure("Color", "w");
for bpi = 1:size( behData.Data, 3 )
    ax = subplot(2,2,bpi);
    if bpi == 2
        auxStack = -squeeze( behData.Data(:,:,bpi) )';
    else
        auxStack = squeeze( behData.Data(:,:,bpi) )';
    end
    auxStack = ( auxStack - median( auxStack, 2) ) ./ max( abs( auxStack ) );
    imagesc( txb*1e3, [],  auxStack ); xline(0, 'k');
    xline( [20, 120], 'LineWidth', 1, 'Color', 'b')
    xlabel('Time [ms]'); ylabel('Trials'); title( behNames(bpi) )
    set( ax, "Box", "off", "Color", "none" )
end
saveFigure(fig, fullfile(figure_path, ...
    "Beh V-0.45 - 0.50 s R25.00 - 350.00 ms", ...
    "All trials all body parts normalised"), true)
clearvars auxStack
cb = colorbar(ax, "Box", "off", "Location", "south");
cb.Ticks = [-1, 1] * 0.85; cb.Label.String = "\leftarrow Direction \rightarrow";
cb.TickLabels = {'Backward', 'Forward'};
linkaxes( findobj(gcf, "Type", "Axes"), "xy")

%% RMS
respWin = sscanf( aInfo.Evoked, "R%f - %f ms")' * 1e-3;
% respWin = [30, 400]*1e-3;
sponWin = -flip(respWin);
[Nb, Nt, Ns] = size( behData.Data );
n = getHesseLineForm([1,0]);

% Spontaneous window
sponFlag = txb > sponWin;
sponFlag = xor( sponFlag(:,1), sponFlag(:,2) );

% Responsive window
evokFlag = txb > respWin;
evokFlag = xor( evokFlag(:,1), evokFlag(:,2) );

respWin_i = [0, diff( respWin )] + 0.12; 
if respWin_i(2) > vwin(2) 
    respWin_i(2) = vwin(2);
end
evokFlag_i = txb > respWin_i;
evokFlag_i = xor( evokFlag_i(:,1), evokFlag_i(:,2) );

evokFlags = [evokFlag, evokFlag_i];

myRMS = @(x) vecnorm(x, 2, 1) ./ size( x, 1 );
funcs = {@(x) x, @(x) diff(x, 1, 1) };
app = [ repmat("", 1,4); repmat( "diff", 1, 4) ];

Nfgs = numel(funcs);
figs = gobjects(Nfgs, 2);
for cf = 1:Nfgs
    figs(cf, 1) = figure( "Color", "w" );
    figs(cf, 2) = figure( "Color", "w" );
    for bpi = 1:size( behData.Data, 3 )

        rwSub = 1;
        ax = subplot( 2, 2, bpi, "Box", "off", ...
            "Color", "none", "Parent", figs(cf, 1) );
        if bpi == 1
            ylabel(ax, 'Evoked')
            rwSub = 2;
        elseif bpi == 4
            xlabel(ax, 'Spontaneous_{RMS}')
        end
        
        aux_x = myRMS( funcs{cf}(behData.Data( sponFlag, :, bpi ) ) );
        aux_y = myRMS( funcs{cf}( ...
            behData.Data( evokFlags(:, rwSub), :, bpi ) ) );
        
        line(ax, aux_x, aux_y, "LineStyle", "none" ); 
        title(ax, join( [behNames(bpi), app(cf,bpi)] ) );
        set( get(ax, "XAxis"), "Scale", "log"); 
        set( get(ax, "YAxis"), "Scale", "log")
        line(ax, xlim(ax), xlim(ax), "LineStyle", "--", ...
            "Color", 0.45*ones(1,3) );
        xticklabels( ax, xticks(ax) ); 
        yticklabels( ax, yticks(ax) );
        text( ax, aux_x, aux_y, num2str( (1:Nt)' ), ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle')
        set(ax, "Box", "off", "Color", "none", "XGrid", "on", ...
            "YGrid", "on")
        Dyx = [aux_x(:), aux_y(:)] * n;
        [~, d_centre, d_scale] = zscore( double( Dyx(pairedStimFlags(:,1)) ) );
        
        ax = subplot(2, 2, bpi, "Box", "off", "Parent", figs(cf, 2) );
        % boxchart(ax, pairedStimFlags * (1:Nccond)', Dyx, "Notch", "on" )
        boxchart( ax, pairedStimFlags * (1:Nccond)', Dyx, ...
            "Notch", "on", "JitterOutliers", "on", "MarkerStyle", "." )
        yyaxis(ax, "right"); line( ax, pairedStimFlags * (1:Nccond)', ...
            my_zscore(Dyx, d_centre, d_scale), "LineStyle", "none")
        yyaxis(ax, "left"); yline(ax, d_centre, 'k' )
        title(ax, behNames(bpi) + app(cf,bpi) );
        set( ax, axOpts{:} ); % set( ax.YAxis, "Scale", "log" )
        xticks(ax, 1:Nccond ); xticklabels( ax, behNames )
        
    end
    saveFigure( figs(cf, 1), fullfile(figure_path, ...
            "Beh V-0.45 - 0.50 s R25.00 - 350.00 ms", ...
            join(["RMS", app(cf,1)]) ), true )

    saveFigure( figs(cf, 2), fullfile(figure_path, ...
            "Beh V-0.45 - 0.50 s R25.00 - 350.00 ms", ...
            join(["RMS", app(cf,2), "boxplots"]) ), true )
end
clearvars aux_* ax figs

%%
sponWeight = (1:sum(sponFlag))/sum(1:sum(sponFlag));

% ln1 = 1:(round( sum( evokFlag )/3 ) );
% ln1 = padarray(ln1, [0, sum( evokFlag ) - numel( ln1 )], "replicate", "post");
ln1 = log10( 1:sum(evokFlag) );

ln2 = sum( evokFlag ):-1:1;
evokWeight = ln1 .* ln2; evokWeight = evokWeight / sum( evokWeight );
rwi = 1;

for cb = 1:size( behData.Data, 3 )
    w_smu = sponWeight*behData.Data(sponFlag,:,cb);
    if cb == 1
        rwi = 2;        
    end
    w_emu = evokWeight*behData.Data( evokFlags(:,rwi), :, cb );
    X = ones( Nt, 1) * txb';
    Y = (1:Nt)' * ones( 1, Nb);
    figure; surf( X, Y, squeeze( behData.Data(:, :, cb) )', ...
        squeeze( behData.Data(:, :, 1) )', "EdgeColor", "interp", ...
        "FaceColor", "none", "EdgeAlpha", 1/3); colormap(gray)
    hold on; line( zeros(Nt,2) + [-0.125,0.03] , 1:Nt, ...
        [w_smu(:), w_emu(:)], ...
        "Marker", "x", "LineStyle", "none", "LineWidth", 2)
    title( behNames(cb) )
    figure; scatter( w_smu, w_emu, [], [w_smu(:), w_emu(:)] * n );
    colormap(-roma+1); colorbar();
    line( xlim, xlim, 'LineStyle', '--', 'Color', 0.45*ones(1,3))
    title( behNames(cb) )
    rwi = 1;
end

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
bxObj = boxchart( reshape( ones(sum( Na ), 1) * (1:4), [], 1), ...
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

%%
video_paths = dir( fullfile( beh_path, "roller*.avi" ) );
vidObj = arrayfun(@(x) VideoReader( expandPath( x ) ), video_paths, fnOpts{:} );

readCSV = @(x) readtable(x, "Delimiter", ",");
fid_paths = dir( fullfile( beh_path, "FrameID*.csv") );
vidTx = arrayfun(@(x) readCSV( expandPath( x ) ), fid_paths, fnOpts{:} );
vidTx = cellfun(@(x) x.Var2/1e9, vidTx, fnOpts{:} ); % nanoseconds

dlcFiles = dir( fullfile( beh_path, "roller*filtered.csv" ) );
fid_paths = dir( fullfile( beh_path, "FrameID*.csv" ) );

exp_path = getParentDir( beh_path, 1);
tf_paths = dir( fullfile( exp_path, "**", "TriggerSignals*.bin") );
fsf_path = dir( fullfile( exp_path, "**", "*_sampling_frequency.mat") );

fs_ephys = load( expandPath( fsf_path ), "fs" ); fs_ephys = fs_ephys.fs;
Ns_intan = [tf_paths.bytes]' ./ 4; % 2 signals x 2 bytes per sample.
Texp_ephys = Ns_intan ./ fs_ephys;

Texp_vid = cellfun(@(x) diff( x([1,end]) ), vidTx );

delta_tiv = Texp_ephys - Texp_vid;

%%
dlcTables = arrayfun(@(x) readDLCData(expandPath(x)), dlcFiles, fnOpts{:} );

