function [t] = createtiles(f, rows, columns)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
t = tiledlayout( f, rows, columns, 'TileSpacing', 'Compact', ...
    'Padding', 'tight');
end