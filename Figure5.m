f = figure("Color", "w"); t = createtiles(f, 2, 1);
tocol = @(x) x(:);
Nb = 8;
clrMap = ember(Nb);
bodypart_names = ["Stim-whisker mean", "Stim-whisker fan arc", ...
    "Nonstim-whisker mean", "Nonstim-whisker fan arc", "Interwhisk arc", ...
    "Symmetry", "Nose", "Roller speed"];
bpn_abb = {'SWM', 'SWF', 'NWM', 'NWF', 'WA', 'S', 'N', 'RS'};
ax = gobjects(2,1);
lsObjs = gobjects( Nb, 2);
xLbls = ["20 - 50 ms", "50 - 200 ms"];
for cx = 1:2
    ax(cx) = nexttile(t);
    cleanAxis( ax(cx) );
    ylabel( ax(cx), xLbls(cx) )
    ytickangle( ax(cx), 90 )
    for cbp = 1:Nb
        lsObjs(cbp,cx) = line(ax(cx), rsq_dt(:,cbp+2), rst_reg_dt(:,3+cx),...
            'LineStyle', 'none', 'Marker', '.', 'Color', clrMap(cbp,:) );
    end
    if cx == 1
        set( get( ax(cx), 'XAxis' ), 'Visible', 'off' )
    end
    legend( ax(cx), lsObjs(:, cx), bpn_abb, "AutoUpdate", "off", ...
        "Box", "off", "Color", "none", "Location", "best", "NumColumns", 2 );
    set( ax(cx), 'TickDir', 'out' )
end
xlabel( ax(cx), "$R^2$", 'Interpreter', 'latex' );
%%
[coeff, score, latent, tsquared, explained] = pca( zscore( rsq_dt(:,3:end), 0, 1 ) );
%%
f = figure( "Color", "w" ); t = createtiles( f, 1, 1 );
ax = nexttile( t ); 
biplot( ax, coeff(:,1:2), "Scores", score(:,1:2), "VarLabels", bpn_abb )
cleanAxis( ax ); title(ax, 'PCA of $R^2$ for all body parts', 'Interpreter', 'latex' )
%% 
f = figure( "Color", "w" ); t = createtiles( f, 2, 1 );
ax = gobjects(2,1); ax(1) = nexttile( t ); hold( ax(1), 'on' );
rsq_dr = (rsq_dt(:,3:end) * coeff(:,1))/2;
inv_dr = (rsq_dt(:,3:end) * coeff(:,2));
line(ax(1), rst_reg_dt(:,4), rsq_dr,...
    'Marker', '.', 'Color', 'k', 'LineStyle', 'none', 'MarkerSize', 16 )
mdl1 = fit_poly( rst_reg_dt(:,4), rsq_dr , 1 );
mdl2 = fit_poly( rst_reg_dt(:,4), rsq_dr , 2 );
prop_ax = 0:0.01:1;
line( ax(1), prop_ax(:), (prop_ax(:).^[2,1,0])*mdl2, 'Color', 'r' );
line( ax(1), prop_ax(:), (prop_ax(:).^[1,0])*mdl1, 'Color', 'b' );
set( get( ax(1), 'XAxis' ), 'Visible', 'off' ); ylabel(ax(1), 'PC1_{R² of all signals}')
ytickangle( ax(1), 90 ); cleanAxis(ax(1));
% mdl2_lm = fitlm( rst_reg_dt(:,4), rsq_dr , "quadratic" );
% coeffDist = getCDFromLM( mdl2_lm );
% CI = createCIfromCD( coeffDist, 0:0.01:1 );
% line( ax(1), (0:0.01:1), CI, 'Color', 'r', 'LineStyle', '--' )
ax(2) = nexttile( t );
line(ax(2), rst_reg_dt(:,4), inv_dr,...
    'Marker', '.', 'Color', 'k', 'LineStyle', 'none', 'MarkerSize', 16 )
ytickangle( ax(2), 90 )
ylabel( ax(2), '_{(SWM,SWF,RS,WA)}\leftarrowPC2\rightarrow_{(N,S,NWM,NWF)}' )
yline( ax(2), 0, 'k--' )
xlabel( ax(2), 'Proportion of responsive units_{20 - 50 ms}')
cleanAxis( ax(2) );
set( ax, 'TickDir', 'out' )
set( f, 'UserData', {coeff, rst_reg_dt, rsq_dt, rsq_dr, inv_dr} )

%%
