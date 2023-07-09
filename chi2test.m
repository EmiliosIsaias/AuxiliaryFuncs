function [p, chiVal] = chi2test(contingencyTbl)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
total = sum(contingencyTbl, "all");
sumCol = sum(contingencyTbl, 1);
sumRow = sum(contingencyTbl, 2);
expectedVal = (sumRow * sumCol)./total;
chiVals = ((contingencyTbl - expectedVal).^2)./expectedVal;
chiVal = sum(chiVals, "all", "omitnan");
df = prod(size(contingencyTbl)-1);
p = chi2cdf(chiVal, df, 'upper');
end