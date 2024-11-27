function coeffDist = getCDFromLM( lmObj, varargin )

p = inputParser;
addRequired(p, 'lmObj', @(x) isa( x, 'LinearModel' ) )
addOptional(p, 'cutoff', 0.05, @(x) isnumeric(x) && numel(x)==1 && x>0 && x<1 )

parse( p, lmObj, varargin{:} )

lmObj = p.Results.lmObj;
cutoff = p.Results.cutoff;

CIBounds = coefCI( lmObj, cutoff );

divFact = fminbnd(@(x) abs( cdf( "Normal", x, 0, 1, "upper" ) - ...
    (cutoff/2) ), -8, 8 );
mu = lmObj.Coefficients.Estimate;
sigma = (CIBounds(:,2) - mu)/divFact;

coeffDist = arrayfun(@(x,y) makedist("Normal", "mu", x, "sigma", y ), mu, sigma);

end