fprintf("Bin size: %.2f\n", configStructure.BinSize_s*1e3)
[psth, txb] = getPSTH_perU_perT(relativeSpkTmsStruct, configStructure);

[Ntrials, Nbins, Nunits] = size( psth{1} );
Nrows = prod( size( psth{1}, [2,3] ) );

neuron = zeros(Nrows, 1, 'uint8');
bin = zeros(Nrows, 1, 'uint16');
counts = zeros(Nrows, Ntrials, 'uint8');

ri = 1;
for cu = 1:size(psth{1}, 3)
    idx = ri:ri+Nbins-1;
    neuron(idx) = cu;
    bin(idx) = 1:Nbins;
    counts(idx,:) = squeeze( psth{1}(:,:, cu) )';
    ri = ri + Nbins;
end