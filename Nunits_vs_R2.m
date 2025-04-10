roller_path = "Z:\Emilio\SuperiorColliculusExperiments\Roller";
expandPath = @(x) fullfile( x.folder, x.name);
bpn_abb = {'SWM', 'SWF', 'NWM', 'NWF', 'WA', 'S', 'N', 'RS'};
fnOpts = {'UniformOutput', false};
scattOpts = {'LineStyle', 'none',  'Color', 'k', 'Marker', '.'};
%%
irn_paths = dir( fullfile( roller_path, "Batch*", "MC", "GADi*", "*", ...
    "ephys*", "Results", "Res *RW20.00-200.00*.mat" ) );

irn_folders = {irn_paths.folder}';
Nm = numel( mice );
Nspm = arrayfun(@(x) numel( x.Sessions ), mice ); % Number of sessions per mouse
Ncl = zeros( sum( Nspm ), 1, 'single' );
Nru = Ncl;
rsp = zeros( size( Ncl, 1 ), 8, 'single' );
rs_vars2load = {'Results', 'Counts'};
cexp = 1;
for cm = 1:Nm
    for cs = 1:Nspm(cm)
        f2load = contains( irn_folders, mice(cm).Name ) & ...
            contains( irn_folders, mice(cm).Sessions(cs).Date ) & ...
            contains( irn_folders, mice(cm).Sessions(cs).Depth );
        load( expandPath(irn_paths(f2load) ), rs_vars2load{:} )
        Ncl(cexp) = size( Results(1).Activity(1).Pvalues, 1 );
        cp_flag = contains( {Results.Combination}', '1 1' );
        Nru(cexp) = nnz( Results(cp_flag).Activity(1).Pvalues < 0.05 );
        rsp(cexp,:) = mice(cm).Sessions(cs).DataTable.R_2_p_L{1}(1,:);
        cexp = cexp + 1;
    end
end
%%
[coeff, score, ~, distTS, explainedVar] = pca( rsp );
ur_mdl = fit_poly( Ncl, score(:,1), 1 );
ur = [0;max(Ncl)].^[1,0] * ur_mdl;
p12_mdl = fit_poly( score(:,1), score(:,2), 1 );

f = figure("Color", "w" ); t = createtiles( f, 3, 2 );
ax(1) = nexttile( t, [2, 2] );
biplot( ax(1), coeff(:,1:2), "Scores", score(:,1:2), "VarLabels", bpn_abb );
grid( ax(1), 'off' ); 
delete( ax(1).Children(2) )
ax(1).YAxisLocation = "origin";
ax(1).YAxis.TickDirection = 'in';
ax(1).XAxisLocation = "origin";
ax(1).XAxis.TickDirection = 'out';
axis( ax(1), 'padded' )

ax(2) = nexttile( t );
line( Ncl, score(:,1), scattOpts{:} ); hold on; 
line( [0; max(Ncl)], ur, 'Color', 'r' )
ylabel(ax(2), 'PC1'); xlabel( ax(2), 'Recorded units' )

ax(3) = nexttile( t );
line( Ncl, score(:,2), scattOpts{:} ); hold on; 
yline( 0, 'Color', 0.15*ones(1,3), 'LineStyle', '--' )
ylabel(ax(3), 'PC2'); xlabel( ax(3), 'Recorded units' )
cleanAxis( ax ); set( ax(2:3), 'TickDir', 'out' )