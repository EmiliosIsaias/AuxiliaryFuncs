function [axOpts] = cleanAxis(ax)
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here
axOpts = {'Box','off','Color','none'};
set( ax, axOpts{:} )
end