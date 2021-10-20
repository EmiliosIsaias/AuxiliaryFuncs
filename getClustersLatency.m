function [outputArg1,outputArg2] = getClustersLatency(spkTms, fs)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%% 
fnOpts = {"UniformOutput", false};
posSpkTms = arrayfun(@(x) cellfun(@(y) y(y > 0), x, fnOpts{:}), spkTms);
fstSpkTms = arrayfun(@(x) cellfun(@(y) y(nnz(y) > 0), x, fnOpts{:}), posSpkTms);
fstPerCls = arrayfun(@(x) cat(2, fstSpkTms{x,:}), (1:size(fstSpkTms,1))',...
    fnOpts{:});

end