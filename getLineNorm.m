function [ n, d ] = getHesseLineForm( mdl )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
angl = angleBetweenLines(0,mdl(1),'rad');
n = [-sin(angl);...
    cos(angl)];
x = 0;
y = x*mdl(1) + mdl(2);
d = [x y]*n;
end

