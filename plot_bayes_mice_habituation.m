histOpts = {'BinLimits', [0, 1], 'NumBins', 32, ...
    'Normalization', 'probability'};
fnOpts = {'UniformOutput', false};
for cm = unique(mice_habituation_cc(:,1))'
    for cs = unique(mice_habituation_cc(mice_habituation_cc(:,1)==cm, 2))'
        cIdx = mice_habituation_cc(:,1)==cm & ...
            mice_habituation_cc(:,2) == cs;
        hats_ms = hats(:, cIdx);
        bH = arrayfun(@(c) histcounts(hats_ms(c,:)*B_scale + B_centre, ...
            histOpts{:}), 1:sum(cIdx), fnOpts{:}); bH = cat(1, bH{:})';
        figure; 
        imagesc(-0.5:0.5:3, 1/64:1/32:1, bH, "Interpolation", "bilinear"); 
        axis('xy'); title(sprintf("Mouse %d, session %d", cm, cs))
        xlabel('Puff intensities'); ylabel('Sessions')
    end
end