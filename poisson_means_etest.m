function [p_value, h] = poisson_means_etest(data1, data2, varargin)
%POISSON_MEANS_ETEST compares two sets of counts assuming they follow a
%Poisson distribution. This test evaluates the null hypothesis that the
%difference between the means of two Poisson distributions is zero.
%   [p_value, h] = poisson_means_etest(data1, data2, varargin)
%% Input parsing
p = inputParser;
% Checking for non-negative integer values
checkData = @(x) all( x>=0 ) && all( mod( x, 1 ) == 0 );
% Checking for alpha
checkAlpha = @(x) isscalar( x ) && isnumeric( x ) && all( x <= 1 & x>0 );

addRequired(p, 'data1', checkData)
addRequired(p, 'data2', checkData)
addOptional(p, 'alpha', 0.05, checkAlpha );

parse(p, data1, data2, varargin{:} )
data1 = p.Results.data1;
data2 = p.Results.data2;
alpha = p.Results.alpha;
%% Perform test

n1 = numel( data1 );
n2 = numel( data2 );
lambda1 = poissfit( data1 );
lambda2 = poissfit( data2 );

z = ( lambda1 - lambda2 ) / sqrt( (lambda1/n1) + (lambda2/n2) );

p_value = 2 * (1 - normcdf( abs( z ) ));
h = p_value < alpha;

end