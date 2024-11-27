function [CI] = createCIfromCD( coeffDist, Xnew, varargin )
%CREATECIFROMCD creates a confidence interval for expected data at alpha
%level
%   CI = createCIfromCD( coeffDist, Xnew );
fnOpts = {'UniformOutput', false};
p = inputParser;
addRequired(p, 'coeffDist');
addRequired(p, 'Xnew', @(x) isnumeric(x) && isvector(x) );
addParameter(p, 'alpha', 0.05, @(x) numel(x) == 1 && isnumeric(x) && x>0 && x<1);

parse(p, coeffDist, Xnew, varargin{:} )
coeffDist = p.Results.coeffDist;
Xnew = p.Results.Xnew;
alpha = p.Results.alpha;

randCoefs = arrayfun(@(x) random( x, [1e4,1] ), coeffDist, fnOpts{:} );
randCoefs = cat(2, randCoefs{:} );
randData = ( Xnew(:).^(0:(numel(coeffDist)-1)) ) * randCoefs';
CI = quantile( randData, (alpha/2)*[1,-1] + [0,1], 2 );

end