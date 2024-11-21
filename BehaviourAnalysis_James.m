%% Analysing behaviour alone

%roller_path = "Z:\Emilio\SuperiorColliculusExperiments\Roller";
% exp_path = ...
%     fullfile( roller_path, "Batch2_ephys/MC/GADi18/211205_C_2450" );
roller_path = fullfile("Z:\James\Practice Data\m7 (training)\");
exp_path = ...
    fullfile( roller_path, "behavior/" );
% Figure overwrite flag
fowFlag = false;
% Annonymus function
expandPath = @(x) fullfile( x.folder, x.name);
% Milli factor
m = 1e-3;
% Looks in the experiment path for objects called ephys* (* is a wild card)
eph_path = dir( fullfile( exp_path, "ephys*" ) );
% Delete everything that is not a folder
eph_path( ~[eph_path.isdir] ) = [];
% Looking folder 'Behaviour'
beh_path = fullfile( exp_path, "Behaviour" );
% If the ephys path is empty, look for *analysis.mat file in the behaviour
% folder.
if ~isempty( eph_path )
    eph_path = expandPath( eph_path );
    figure_path = fullfile( eph_path, "Figures" );
    af_path = dir( fullfile( eph_path , "*analysis.mat" ) );
    af_path = expandPath( af_path );
elseif exist( beh_path, "dir" )
    af_path = expandPath( dir( fullfile( beh_path, "*analysis.mat") ) );
    figure_path = fullfile( beh_path, "Figures" );
else
    % If there is no organisation in the folders, look for the file in the
    % 'root' folder a.k.a. experiment folder. Then give up.
    beh_path = exp_path;
    af_path = expandPath( dir( fullfile( beh_path, "*analysis.mat") ) );
    figure_path = fullfile( beh_path, "Figures" );
end

[~, af_name] = fileparts( af_path );
expName = extractBefore(af_name, "analysis");
load( af_path, "Conditions", "fs")

fnOpts = {'UniformOutput', false};
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
consCond = 1;
Nccond = length( consCond );
%prmSubs = nchoosek( 1:Nccond, 2 );

% pairedStimFlags = arrayfun(@(c) any( ...
%     Conditions(1).Triggers(:,1) == ...
%     reshape( Conditions(c).Triggers(:,1), 1, [] ), 2 ), consCond, fnOpts{:} );
% pairedStimFlags = cat( 2, pairedStimFlags{:} );

consCondNames = string( { Conditions( consCond ).name  } );

% [behRes, behFig_path, behData, aInfo] = analyseBehaviour( beh_path, ...
%     "ConditionsNames", cellstr( consCondNames ), ...
%     "PairedFlags", pairedStimFlags, ...
%     "FigureDirectory", figure_path, ...
%     "ResponseWindow", [25, 350] * m, ...
%     "ViewingWindow", [-450, 500] * m, ...
%     "figOverWrite", fowFlag );
[behRes, behFig_path, behData, aInfo] = analyseBehaviour( beh_path, ...
    "ConditionsNames", cellstr( consCondNames ), ...
    "FigureDirectory", figure_path, ...
    "ResponseWindow", [25, 350] * m, ...
    "ViewingWindow", [-450, 500] * m, ...
    "figOverWrite", fowFlag );

if ldFlag
    load( expandPath( dir( fullfile( beh_path, "RollerSpeed*.mat" ) ) ), "fr")
end
[Ns, Nt, Nb] = size( behData.Data );

vwin = sscanf( aInfo.VieWin, "V%f - %f s")';
mdlt = fit_poly( [1, Ns], vwin + [1,-1] * (1/(2 * fr) ), 1 );
txb = ( (1:Ns)'.^[1,0] ) * mdlt;
behNames = string( { behRes(1).Results.BehSigName } );

%% Normalised amplitud by absolute maximum

rollYL = "Roller speed [cm/s]";
yLabels = [repmat("Angle [°]", 1, Nb-1), rollYL];
sym_flag = contains( behNames, "symmetry", "IgnoreCase", true );
yLabels(sym_flag) = "Symmetry [a.u.]";
cbLabels = strings(Nb, 2);
cbLabels([1,3],:) = repmat(["Retract", "Protract"],2,1);
cbLabels([2,4,5],:) = repmat(["Closed", "Opened"],3,1);
cbLabels(6,:) = ["Away", "Puff"];
cbLabels(7,:) = ["Puff", "Away"];
cbLabels(8,:) = ["Backward", "Forward"];
screen_size = get(0, 'ScreenSize' );
pxHeight = screen_size(4)*0.7;
lowBound = screen_size(4)*1/5;

possCols = [2,3,5];
Ncols = possCols( find( mod( Nb, possCols ) == 0, 1, 'first' ) );
Nrows = Nb / Ncols;

newAx = @(f) subplot( Nrows, Ncols, ix, "NextPlot", "add", "Parent", f );

createtiles = @(f) tiledlayout( f, Nrows, Ncols, ...
    'TileSpacing', 'Compact', 'Padding', 'tight');

isendrow = @(ix) ( (ix/Ncols) + 1) > Nrows;

newFigure = @(x,y,fn) figure( "Color", "w", "Position", ...
    [x, y, pxHeight/sqrt(2), pxHeight], "Name", fn );

fig = newFigure(0, lowBound, "");

t = createtiles( fig );
for cbi = 1:Nb
    ax = nexttile(t);
    auxStack = squeeze( behData.Data(:,:,cbi) )';
    if cbi ~= 6
        auxStack = ( auxStack - median( auxStack, 2) );
        auxStack = auxStack ./ max( abs( auxStack ), [], 2 );
    end
    imagesc( ax, txb*1e3, [],  auxStack ); xline(ax, 0, 'k');
    xline( ax, [20, 120], 'LineWidth', 1, 'Color', 'b')
    title( ax, behNames(cbi) )
    if mod( cbi, Ncols ) == 1
        ylabel(ax, 'Trials'); 
    else
        ax.YAxis.Visible = 'off';
    end
    if isendrow( cbi )
        xlabel(ax, 'Time [ms]'); 
    else
        ax.XAxis.Visible = 'off';
    end
    set( ax, axOpts{:} )
    
    cb = colorbar(ax, "Box", "off", "Location", "west");
    cb.Ticks = [min(auxStack(:)), max(auxStack(:))] * 0.85;
    cb.Label.String = yLabels(cbi);
    cb.TickLabels = cellstr(cbLabels(cbi,:));
    cb.TickDirection = "none";
end
axis( findobj( fig, "Type", "Axes" ), ...
    [1e3*txb([1,end])', [1,Nt] + [-1,1]*(1/2)] )
%cb.TickLabels = {'Backward', 'Forward'};
linkaxes( findobj(fig, "Type", "Axes"), "xy")

saveFigure(fig, fullfile( behFig_path, ...
    "All trials all body parts normalised" ), true, fowFlag)

clearvars auxStack

%% L-norms and derivative
my_zscore = @(x, m, s) ( x - m ) ./ ( s .* (s~=0) + 1 .* (s==0) );

respWin = sscanf( aInfo.Evoked, "R%f - %f ms")' * 1e-3;
% respWin = [30, 400]*1e-3;
sponWin = -flip(respWin);
n = getHesseLineForm([1,0]);

% Spontaneous window
sponFlag = txb > sponWin;
sponFlag = xor( sponFlag(:,1), sponFlag(:,2) );

% Responsive window
respWin_aux = respWin;
if respWin(1) < 0.12
    respWin_aux = respWin + 0.1;
end

if respWin_aux(2) > vwin(2)
    respWin_aux(2) = vwin(2);
end
respWin = [repmat(respWin_aux,2,1); repmat( respWin, Nb-2, 1)];

evokFlags = arrayfun(@(x) txb > respWin(x,:), 1:Nb, fnOpts{:} );
evokFlags = cellfun(@(x) xor(x(:,1), x(:,2) ), evokFlags, fnOpts{:} );
evokFlags = cat( 2, evokFlags{:} );

myNorm = @(x, l) vecnorm(x, l, 1);
funcs = {@(x) x, @(x) diff(x, 1, 1) };
app = [ repmat("", 1,Nb); repmat( "diff", 1, Nb) ];

Nfgs = numel(funcs);
figs = gobjects(Nfgs, 2);

for l = [1, 2, inf]

    for cf = 1:Nfgs
        figs(cf, 1) = newFigure(0, lowBound, "");
        %figure( "Color", "w", "Position", ...
         %   [0, 20, pxHeight/sqrt(2), pxHeight] );
        figs(cf, 2) = newFigure(pxHeight/sqrt(2), lowBound, "");... figure( "Color", "w", "Position", ...
            %[pxHeight/sqrt(2), 20, pxHeight/sqrt(2), pxHeight] );
        t1 = createtiles( figs(cf, 1) );
        t2 = createtiles( figs(cf, 2) );
        for cb = 1:Nb

            ax = nexttile(t1);
            if cb == 1
                ylabel(ax, 'Evoked')
            elseif cb == Nb
                xlabel(ax, 'Spontaneous')
                lgObj = legend( ax, scObj, consCondNames, lgOpts{:}, ...
                "AutoUpdate", "off" );
            end

            aux_x = myNorm( funcs{cf}(behData.Data( sponFlag, :, cb ) ), l );
            aux_y = myNorm( funcs{cf}( ...
                behData.Data( evokFlags(:, cb), :, cb ) ), l );

            scObj = arrayfun(@(c) line(ax, aux_x(pairedStimFlags(:,c)), ...
                aux_y(pairedStimFlags(:,c)), "LineStyle", "none", ...
                "Marker", "." ), 1:Nccond);
            
            title(ax, join( [sprintf("L%d", l), behNames(cb), ...
                app(cf,cb)] ) );
            set( get(ax, "XAxis"), "Scale", "log");
            set( get(ax, "YAxis"), "Scale", "log")
            line(ax, xlim(ax), xlim(ax), "LineStyle", "--", ...
                "Color", 0.45*ones(1,3) );
            xticklabels( ax, xticks(ax) );
            yticklabels( ax, yticks(ax) );
            % text( ax, aux_x, aux_y, num2str( (1:Nt)' ), ...
            %     'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle')
            set(ax, axOpts{:}, "XGrid", "on", "YGrid", "on")
            Dyx = [aux_x(:), aux_y(:)] * n;
            % [~, d_centre, d_scale] = zscore( double( Dyx(pairedStimFlags(:,1)) ) );

            ax = nexttile(t2);
            % boxchart(ax, pairedStimFlags * (1:Nccond)', Dyx, "Notch", "on" )
            trFlag = any( pairedStimFlags, 2 );
            boxchart( ax, pairedStimFlags(trFlag,:) * (1:Nccond)', ...
                Dyx(trFlag), "Notch", "on", "JitterOutliers", "on", ...
                "MarkerStyle", "." )
            % yyaxis(ax, "right"); line( ax, pairedStimFlags * (1:Nccond)', ...
            %     my_zscore(Dyx, d_centre, d_scale), "LineStyle", "none")
            % yyaxis(ax, "left"); 
            yline( ax, 0, 'k' )
            title(ax, join( [sprintf("L%d", l), ...
                behNames(cb), app(cf,cb)] ) );
            set( ax, axOpts{:} ); % set( ax.YAxis, "Scale", "log" )
            xticks(ax, 1:Nccond ); 
            if isendrow( cb )
                xticklabels( ax, consCondNames )
            else
                ax.XAxis.Visible = 'off';
            end

        end
        saveFigure( figs(cf, 1), fullfile(behFig_path, ...
            join([sprintf("L%d norm", l), app(cf,1)]) ), true, fowFlag )

        saveFigure( figs(cf, 2), fullfile(behFig_path, ...
            regexprep( join( [sprintf( "L%d norm", l ), app(cf,2), ...
            "boxplots"] ), ' +', ' ' ) ), true, fowFlag )
    end
    clearvars aux_* ax figs
end

%% Weigthed mean
sponWeight = (1:sum(sponFlag))/sum(1:sum(sponFlag));
ln1 = log10( 1:sum(evokFlags(:,1)) );
ln2 = sum( evokFlags(:,1) ):-1:1;
evokWeight = ln1 .* ln2; evokWeight = evokWeight / sum( evokWeight );

clrMap = lines(Nccond);

xPos = (0:2)*pxHeight/sqrt(2);
Na = sum( pairedStimFlags );
figs = gobjects( 3, 1 );
ts = figs;
figNames = ["Weighted mean"; "Line distance boxplots"; 
    "Vector magnitude boxplots"];

for cf = 1:numel(figs)
    figs(cf) = newFigure( xPos(cf), lowBound, figNames(cf) );
    ts(cf) = createtiles( figs( cf ) );
end


for cb = 1:Nb
    w_smu = reshape( sponWeight * behData.Data(sponFlag,:,cb), [], 1 );
    
    w_emu = reshape( evokWeight * behData.Data( ...
        evokFlags(:,cb), :, cb ), [], 1 );

    
    ax = nexttile( ts(1) ); set(ax, 'NextPlot', 'add' );
    scObj = arrayfun(@(c) scatter(ax, w_smu(pairedStimFlags(:,c)), ...
        w_emu(pairedStimFlags(:,c)), '.', "MarkerEdgeColor", clrMap(c,:) ), ...
        1:Nccond);

    if cb == Nb
        xlabel( ax, 'Spontaneous', 'FontSize', 8 );
        legend( ax, scObj, consCondNames, lgOpts{:}, "AutoUpdate", "off");
    elseif cb == 1
        ylabel( ax, 'Evoked', 'FontSize', 8 );
    end

    title(ax, behNames(cb) );
    set( ax, axOpts{:}, "XAxisLocation", "origin", ...
        "YAxisLocation", "origin" ); grid( ax, "on" ); %axis( ax, 'square' )
    line(ax, xlim(ax), xlim(ax), 'LineStyle', '--', 'Color', 0.45*ones(1,3))

    ax = nexttile( ts(2) ); set(ax, 'NextPlot', 'add' );

    bxObj = boxchart(ax, pairedStimFlags(trFlag,:) * (1:Nccond)', ...
        [w_smu(trFlag), w_emu(trFlag)] * n, ...
        "Notch", "on", "JitterOutlier", "on", "MarkerStyle", ".", ...
        "Boxfacecolor", "k", "Markercolor", "k" );

    if cb == Nb
        legend( ax, bxObj, '$\vec{x} \cdot n + d$', ...
            'Interpreter' ,'latex' , lgOpts{:}, "AutoUpdate", "off" );
    elseif cb == 1
        ylabel( ax, 'Distance from line' )
    end

    xticks( ax, 1:Nccond ); 
    if isendrow(cb)
        xticklabels( ax, consCondNames )
    else
        ax.XAxis.Visible = 'off';
    end
    
    title(ax, behNames(cb) );
    set( ax, axOpts{:} ); yline( ax, 0, 'Color', 0.75*ones(1,3), ...
        'LineWidth', 1/3);

    ax = nexttile( ts(3) ); set(ax, 'NextPlot', 'add' );
    bxObj = boxchart(ax, pairedStimFlags(trFlag,:) * (1:Nccond)', ...
        vecnorm( [w_smu(trFlag), w_emu(trFlag)], 2, 2 ), ...
        "Notch", "on", "JitterOutlier", "on", "MarkerStyle", "." );
    if cb == Nb
        legend( ax, bxObj, 'L-2 norm', ...
            lgOpts{:}, "AutoUpdate", "off", "interpreter", "latex" );
    elseif cb == 1
        ylabel( ax, 'Distance from origin' )
    end

    xticks( ax, 1:Nccond ); 
    if isendrow(cb)
        xticklabels( ax, consCondNames )
    else
        ax.XAxis.Visible = 'off';
    end

    title(ax, behNames(cb) ); %axis( ax, 'square' )
    set( ax, axOpts{:} );
    
    for cc = 1:Nccond
        behRes(cc).Results(cb).Puff_Effect = ...
            [w_smu( pairedStimFlags(:,cc) ), ...
            w_emu( pairedStimFlags(:,cc) )] * n;
        behRes(cc).Results(cb).Baseline_L2 = ...
            vecnorm( [w_smu( pairedStimFlags(:,cc) ), ...
            w_emu( pairedStimFlags(:,cc) )], 2, 2);
    end
end

wmFigNames = ["Weighted mean scatter"; ...
    "Line distance boxplots"; ...
    "Origin distance"];

arrayfun(@(x, f) saveFigure( x, fullfile( behFig_path, f ), true, fowFlag), ...
    figs(:), wmFigNames(:))

clearvars ax figs

%% Amplitude index and trial proportion

if ~exist( "fr", "var" ) && ldFlag
    load( expandPath( dir( fullfile( beh_path, "RollerSpeed*.mat" ) ) ), "fr")
    ldFlag = false;
end

[pAreas, ~, behAreaFig] = createBehaviourIndex(behRes);
behMeasures = string({behAreaFig.Name});
biFigPttrn = behMeasures+"%s";
biFigPttrn = arrayfun(@(s) sprintf(s, sprintf(" %s (%%.3f)", ...
    consCondNames ) ), biFigPttrn );

for it = 1:numel(behMeasures)
    behRes = arrayfun(@(bs, ba) setfield( bs, ...
        strrep( behMeasures(it), " ", "_" ), ba), behRes(:), pAreas(:,it) );
end

arrayfun(@(f) set( f, 'UserData', behRes ), behAreaFig );

biFN = arrayfun(@(s) sprintf( biFigPttrn(s), pAreas(:,s) ), 1:numel(behMeasures) );

arrayfun(@(f, fn) saveFigure(f, fullfile(behFig_path, fn), true, fowFlag), ...
    behAreaFig(:), biFN(:) );



%% Count figure
trMvFlag = arrayfun(@(cr) behRes(1).Results(cr).MovStrucure.MovmentFlags, ...
    1:size(behRes(1).Results,2), fnOpts{:}); trMvFlag = cat(3, trMvFlag{:});
BIscaleMat = sum(trMvFlag,3);
BIscale = arrayfun(@(cc) BIscaleMat(pairedStimFlags(:,cc), cc), 1:Nccond, ...
    fnOpts{:});
hstOpts = {'BinMethod', 'integers', 'BinLimits', [0,Nb] + [-1,1]/2};
[hg, hg_bin] = cellfun(@(c) histcounts(c, hstOpts{:}), ...
    BIscale, fnOpts{:});
hg = cat(1, hg{:}); hg_bin = cat(1, hg_bin{:});

% [p_amp, h_amp] = ranksum(cat(1, zamp{1,:}), cat(1, zamp{2,:}));

clrMap = lines(Nccond);
countFig = figure; ax(1) = subplot(10,1,1:8);
bar(ax(1), (0:Nb)', (hg./sum(hg,2))', 'EdgeColor', 'none'); hold on;
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

saveFigure(countFig, fullfile(behFig_path, countFigName), true, fowFlag);