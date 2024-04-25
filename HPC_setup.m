%%
home_path = '/gpfs/bwfor/home/hd/hd_hd/hd_bf154/';
repo_paths = cellfun(@(x) char( fullfile( home_path, x) ), ...
    {'NeuroNetzAnalysis', 'AuxiliaryFuncs', 'Scripts'},  fnOpts{:} );
addpath( repo_paths{:} )