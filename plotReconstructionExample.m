phOpts = {'EdgeColor', 'none', 'FaceAlpha', 0.25, 'FaceColor'};
m = 1e-3; k = 1e3;
f = figure('Color', 'w'); t = createtiles( f, 4, 2 );
ax = gobjects( 8, 1 );
yLbls = [repmat({'Angle [°]'}, 5, 1 ); {'a.u.'}; {'Angle [°]'};
    {'Speed [\frac{cm}{s}]'}];
for cb = 1:Nb
    ax(cb) = nexttile(t); hold( ax(cb), 'on' );
    title( ax(cb), bodypart_names(cb) )
    line( ax(cb), k*trial_tx, oMu(:,cb), 'Color', clrMap(2,:), 'LineWidth', 2 )
    line( ax(cb), k*trial_tx, rMu(:,cb), 'Color', clrMap(1,:), 'LineWidth', 2 )
    patch( ax(cb), k*trial_tx([1:end,end:-1:1]), ...
        mat2ptch(getSEM(obStack{1}(:,:,cb))), 1, phOpts{:}, clrMap(2,:) )
    patch( ax(cb), k*trial_tx([1:end,end:-1:1]), ...
        mat2ptch(getSEM(rbStack{1}(:,:,cb))), 1, phOpts{:}, clrMap(1,:) )
    if any( cb == [1,3,7])
        yticklabels( ax(cb), yticks( ax(cb) ) - 90 )
    end
    if ~any( cb == [7,8] )
        set( get( ax(cb), 'XAxis' ), 'Visible', 'off' )
    else
        xlabel( ax(cb), 'Time [ms]' )
    end
    if cb == 8
        ylabel( ax(cb), ['$',yLbls{cb},'$'], 'Interpreter', 'latex' )
    else
        ylabel( ax(cb), yLbls{cb} )
    end
    ytickangle( ax(cb), 90 );
end
cleanAxis(ax);
set( ax, 'TickDir', 'out' )
legend( ax(cb), {'Original', 'Reconstructed'}, 'Box', 'off', ...
    'Color', 'none', 'Location', 'best', 'AutoUpdate', 'off' )
linkaxes( ax, 'x'); xlim( ax, [-300, 800] )
set( f, 'UserData', {rbStack, obStack} )