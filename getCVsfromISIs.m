function [CV2, CVsqr] = getCVsfromISIs(ISI)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
CVsqr = var(ISI)./(mean(ISI).^2);
isiVar = abs(diff(ISI));
isiSum = zeros(numel(isiVar),1);
for cisi = 1:numel(isiVar)
    isiSum(cisi) = ISI(cisi+1) + ISI(cisi);
end
CV2 = mean((2.*isiVar)./(isiSum));
end

