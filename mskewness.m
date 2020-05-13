function s = mskewness(x,flag,dim,mns)
if nargin < 2 || isempty(flag)
    flag = 1;
end
if nargin < 3 || isempty(dim)
    % The output size for [] is a special case, handle it here.
    if isequal(x,[]), s = NaN('like',x); return; end;

    % Figure out which dimension nanmean will work along.
    dim = find(size(x) ~= 1, 1);
    if isempty(dim), dim = 1; end
end

% Need to tile the output of nanmean to center X.
tile = ones(1,max(ndims(x),dim));
tile(dim) = size(x,dim);

% Center X, compute its third and second moments, and compute the
% uncorrected skewness.
x0 = x - repmat(mns, tile);
s2 = nanmean(x0.^2,dim); % this is the biased variance estimator
m3 = nanmean(x0.^3,dim);
s = m3 ./ s2.^(1.5);

% Bias correct the skewness.
if flag == 0
    n = sum(~isnan(x),dim);
    n(n<3) = NaN; % bias correction is not defined for n < 3.
    s = s .* sqrt((n-1)./n) .* n./(n-2);
end

end