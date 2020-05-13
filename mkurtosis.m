function k = mkurtosis(x,flag,dim,mns)
if nargin < 2 || isempty(flag)
    flag = 1;
end
if nargin < 3 || isempty(dim)
    % The output size for [] is a special case, handle it here.
    if isequal(x,[]), k = NaN('like',x); return; end;

    % Figure out which dimension nanmean will work along.
    dim = find(size(x) ~= 1, 1);
    if isempty(dim), dim = 1; end
end

% Need to tile the output of nanmean to center X.
tile = ones(1,max(ndims(x),dim));
tile(dim) = size(x,dim);

% Center X, compute its fourth and second moments, and compute the
% uncorrected kurtosis.
x0 = x - repmat(mns, tile);
s2 = nanmean(x0.^2,dim); % this is the biased variance estimator
m4 = nanmean(x0.^4,dim);
k = m4 ./ s2.^2;

% Bias correct the kurtosis.
if flag == 0
    n = sum(~isnan(x),dim);
    n(n<4) = NaN; % bias correction is not defined for n < 4.
    k = ((n+1).*k - 3.*(n-1)) .* (n-1)./((n-2).*(n-3)) + 3;
end
end