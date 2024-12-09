function [rtm] = getClusterRange(coeff, MaxComp, MaxClust )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
y = pdist( coeff(:,1:MaxComp) );
z = linkage( y );
rtm = cluster( z, "maxclust", MaxClust );
end