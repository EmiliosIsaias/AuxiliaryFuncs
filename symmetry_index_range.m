symind = @(x,y) cosd(x) - cosd(y);
[X,Y] = meshgrid( 0:0.5:180 );
Z = symind( X, Y );
f = figure("Color", "w"); t = createtiles(f, 1, 1);
ax = nexttile( t ); contour(ax, X, Y, Z, -2:2 );
cleanAxis( ax ); axis(ax, 'equal' )
xlabel('Stimulated whisker mean [°]')
ylabel('Non-stimulated whisker mean [°]'); ytickangle( ax, 90 )
colormap( [1, 0.4, 0;0,0,0;1, 0.4, 0] )
hold(ax, "on" )
set( ax, 'TickDir', 'out' )