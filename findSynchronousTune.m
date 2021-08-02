function [complexityHist] = findSynchronousTune(relSpkTms, binSzs, tmWin, varargin)
%FINDSYNCHRONOUSTUNES bins the given relative spike times per each bin size
%and creates a histogram for how many simultaneous spikes accour in a given
%bin.
%   Detailed explanation goes here

%% Parse inputs

p = inputParser;
% Required arguments
checkSpkStruct = @(x) any([isstruct(x), isfield(x,{'name','SpikeTimes'})]);

checkBnSzs = @(x) any([isnumeric(x), all(x>0), isrow(x) | iscolumn(x)]);

checkTmWin = @(x) all([numel(x) == 2, x(1) < x(2)]);

p.addRequired('relSpkTms', checkSpkStruct);
p.addRequired('binSzs', checkBnSzs);
p.addRequired('tmWin', checkTmWin);


p.parse(relSpkTms, binSzs, tmWin, varargin{:});

relSpkTms = p.Results.relSpkTms;
binSzs = p.Results.binSzs;
tmWin = p.Results.tmWin;

%% Main loop
fnOpts = {'UniformOutput', 0};
Ncond = numel(relSpkTms);
Ncl = size(relSpkTms(1).SpikeTimes,1);
cplxOpts = {'BinMethod', 'integers', 'BinLimits', [1, Ncl],...
    'Normalization', 'probability'};
Nbns = size(binSzs(:), 1);


spkTms_ClandTr = arrayfun(@(x) arrayfun(@(y) cat(2, x.SpikeTimes{:,y}),...
    (1:size(x.SpikeTimes,2))', fnOpts{:}), relSpkTms(:), fnOpts{:});

complexityHist = zeros(Ncl, Nbns, Ncond, 'single');
for cbn = 1:Nbns
    hstOpts = {'BinWidth', binSzs(cbn), 'BinLimits', tmWin};
    hbn = cellfun(@(x) cellfun(@(y) histcounts(y, hstOpts{:}),...
        x, fnOpts{:}), spkTms_ClandTr, fnOpts{:});
    %hbn = cellfun(@(x) mean(cat(1, x{:})), hbn, fnOpts{:});
    hbn = cellfun(@(x) cat(1, x{:}), hbn, fnOpts{:});
    auxComplex = cellfun(@(x) histcounts(x, cplxOpts{:})', hbn, fnOpts{:});
    complexityHist(:, cbn, :) = cat(3, auxComplex{:});
end

end