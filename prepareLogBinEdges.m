function [binCenters, binEdges, lData, ts] = prepareLogBinEdges(data, nBins)
lData = log10(data); mnVal = min(lData); mxVal = max(lData);
dRange = range(lData); ts = dRange / (nBins - 1);
binCenters = mnVal:ts:mxVal; binDomain = [mnVal;mxVal] + [-1;1]*(ts/2);
binEdges = binDomain(1):ts:binDomain(2);
end

