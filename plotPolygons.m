function [pchObj, ax] = plotPolygons(poly_coords, f, varargin) %#ok<*INUSD>
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here
%% Input parser
checkPC = @(x) isnumeric(x) & ~isreal(x);
checkF = @(x) isa(x,'matlab.ui.Figure');
checkCM = @(x) ismatrix(x) & size( x, 2 ) == 3 & size( x, 1 ) >= 1 & ...
    all( x >= 0 & x <= 1, "all" );
p = inputParser;
addRequired(p, 'poly_coords', checkPC )
addRequired(p, 'f', checkF )
addParameter(p, 'clrMap', lines( size( poly_coords, 2 ) ), checkCM )

parse(p, poly_coords, f, varargin{:} )

poly_coords = p.Results.poly_coords;
f = p.Results.f;
clrMap = p.Results.clrMap;
if size( clrMap, 1 ) ~= size( poly_coords, 2 )
    fprintf(1, 'The colormap and the number of polygons given do not match\n')
    fprintf(1, 'colormap: %d , polygons: %d', size( clrMap, 1 ), ...
        size( poly_coords, 2 ) )
    clrMap = lines( size( poly_coords, 2 ) );
end
%% Function
ax = findobj( f, 'Type', 'axes' );

pchObj = arrayfun(@(c) patch( ax, real( poly_coords(:,c) )', ...
    imag( poly_coords(:,c) )', clrMap(c,:), 'EdgeColor', 'none', ...
    'FaceAlpha', 0.5 ), 1:size( poly_coords, 2 ) );
end