function [f, z_axis, poly_coords] = createPolarPlotPolygons( ai )
%CREATEPOLARPLOTPOLYGONS returns figure object, N complex vector axes, and
%the coordinates for the given matrix ai. 
%   [f, z_axis, poly_coords] = createPolarPlotPolygons( ai )

figOpts = checkSystem4Figures();
% createtiles = @(f,r,c) tiledlayout( f, r, c, ...
%     'TileSpacing', 'Compact', 'Padding', 'tight');
[Nb, ~] = size( ai );

radAxis = (0:Nb-1)*(2*pi/Nb);
z_axis = exp(1i*radAxis(:));

poly_coords = ai .* z_axis;

f = figure('Color', 'w', figOpts{:} );
t = createtiles( f, 1, 1 ); ax = nexttile(t);

% Drawing polar axis 
zend = Nb/2;
if mod( Nb, 2 )
    % Odd number of axes
    zend = Nb;
end
line(ax, [-real(z_axis(1:zend)), real(z_axis(1:zend))]', ...
    [-imag(z_axis(1:zend)), imag(z_axis(1:zend))]', 'LineWidth', 0.1, ...
    'Color', 0.45*ones(1,3));
arrayfun(@(x) rectangle(ax, 'Position', [repmat(-x,1,2), repmat(x*2,1,2)], ...
    'Curvature', [1,1], 'EdgeColor', 0.45*ones( 1, 3 ) ), 1:-0.25:0.25 );
text( ax, 0.25:0.25:1, zeros(1,4), string( ( 1:4 )'/4 ), ...
    "HorizontalAlignment", "left", "VerticalAlignment", "cap" )
cleanAxis(ax); axis(ax, 'off', 'equal');
end