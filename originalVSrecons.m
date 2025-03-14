%% Comparison observed vs reconstructed
% Assuming reconstructed and observed already in workspace.

[~, Nbs, Nexp] = size( r2_res_c );
getSEM = @(x, d) cat( d, mean( x, d, "omitmissing" ), ...
    std( x, 1, d ,"omitmissing" )./sqrt( size( x, d ) ) );
mat2ptch = @(x) [x(1:end,:)*[1;1]; x(end:-1:1,:)*[1;-1]];
bp_names = ["Stim-whisker mean", "Stim-whisker fan arc", ...
    "Nonstim-whisker mean", "Nonstim-whisker fan arc", ...
    "Interwhisker arc", "Symmetry", "Nose", "Roller speed"];
yLbls = [repmat({'Angle [\circ]'}, 5, 1 ); {'a.u.'}; {'Angle [\circ]'};
{'Speed [$\frac{cm}{s}$]'}];
trial_tmdl = fit_poly( [1,params.Nb], ...
    params.relative_window + [1,-1]*(params.bin_size/2), 1 );
trial_tx = ( (1:params.Nb)'.^[1,0] ) * trial_tmdl;

f = figure("Color", "w"); t = createtiles( f, Nbs/2, 2 );
ax = gobjects( Nbs, 1 );
laser_blue = [0, 132, 255]/255;
recons_green = [32, 138, 0]/255;
clrMap = [0.15*ones(1,3); recons_green];

behSignals = {y_trials, y_ptrials};
extras = {'Interpreter', 'tex'};
for cbs = 1:Nbs
    ax(cbs) = nexttile( t ); hold( ax(cbs), "on" );
    for cs = 1:2        patch( ax(cbs), trial_tx([1:end,end:-1:1]), mat2ptch( ...
            getSEM( behSignals{cs}(:,:,cbs), 2 ) ), 1, 'EdgeColor', 'none', ...
            'facecolor', clrMap(cs,:), 'FaceAlpha', 1/2 )
        line( ax(cbs), trial_tx, mean( behSignals{cs}(:,:,cbs), 2 ), ...
            'Color', clrMap(cs,:), 'LineWidth', 2 );
    end
    title( ax(cbs), bp_names(cbs) )
    set( ax(cbs), 'TickDir', 'out' );
    ytickangle( ax(cbs), 90 );
    if cbs == Nbs
        extras{2} = 'latex';
    end
    ylabel( ax(cbs), yLbls{cbs}, extras{:} )
end
axOpts = cleanAxis( ax );
arrayfun(@(x) set( get( x, "XAxis" ), "Visible", "off" ), ax(1:6) )
linkaxes( ax, 'x' )
arrayfun(@(x) yticklabels( x, yticks(x) - 90 ), ax([1,3,7]) )
xlim( ax, [-0.4, 0.8] );
arrayfun(@(x) xticklabels( x, round(xticks(x),1)*1e3 ), ax(7:8) )
xlabel( ax(end), 'Time [ms]' )

axCh = get( ax(1), 'Children' );
axCh(~arrayfun(@(x) isa(x, 'matlab.graphics.primitive.Line'), axCh )) = [];
legend( ax(1), axCh, {'Reconstructed','Observed'}, axOpts{:}, ...
    'Location', 'best', 'Autoupdate', 'off' )
%%
arrayfun(@(x) patch( x, m*[-100,200,200,-100], ...
    tocol( repmat( (ylim(x) * [1;0]) * [1,0.98], 2, 1 ) ), ...
    [0, 132, 255]/255, "EdgeColor", "none", "FaceAlpha", 0.5 ), ax )