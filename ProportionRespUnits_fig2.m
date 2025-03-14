% Responsive units in all experiments between 20 and 200 ms
roller_path = "Z:\Emilio\SuperiorColliculusExperiments\Roller";
res_paths = dir( fullfile( roller_path, "Batch*", "*", "*", "*", "ephys*", ...
    "Results", "Res VW* ms RW20.00-200.00 ms* PuffAll.mat" ) );
map_paths = dir( fullfile( roller_path, "Batch*", "*", "*", "*", ...
    "ephys*", "Results", "Map *.mat" ) );
% {map_paths( ~ismember({map_paths.folder}', {res_paths.folder}' ) ).folder}'
for 1:numel( res_paths )
    
end