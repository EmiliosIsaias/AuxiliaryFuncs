function [results, f] = AnBeh_Bypass(data_path, resWin, varargin)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%% Input validation

checkDir = @(x) (ischar(x) || isstring(x)) && exist(x, 'dir');
checkRW = @(x) isnumeric(x) & numel(x) == 2 & (x(2) > x(1));

p = inputParser;

addRequired(p, 'data_path', checkDir)
addRequired(p, 'resWin', checkRW)
addParameter(p, 'sponWin', -flip(resWin), checkRW)

parse(p, data_path, resWin);

data_path = p.Results.data_path;
brWin = p.Results.resWin;
bsWin = p.Results.sponWin;
%% Auxiliary variables and functions
fnOpts = {'UniformOutput', false};

lgOpts = {'Box', 'off', 'Color', 'none', 'Location', 'best'};
vaxOpts = cellstr( ["HorizontalAlignment", "center", ...
    "VerticalAlignment", "baseline", "Rotation"] );
%TODO: Get bodypart names from somewhere else!
bodypart_names = ["Stim-whisker mean", "Stim-whisker fan arc", ...
"Nonstim-whisker mean", "Nonstim-whisker fan arc", "Interwhisk arc", ...
"Symmetry", "Nose", "Roller speed"];

expandName = @(x) fullfile( x.folder, x.name );
m = 1e-3; k = 1e3;
%% Load regression file
regFile = dir( fullfile( data_path, "Regression CW*.mat" ) );
if ~isempty( regFile ) && numel( regFile ) == 1
    load( expandName( regFile ), 'DX', 'mdlAll_ind', 'params')
else
    %TODO: Being able to select a regression for specific parameters
    fprintf(1, 'Regression either didn''t work or have different files\n')
    return
end
% Colormap: grey for original and PCB-green for reconstructed
clrMap = flip([0.15*ones(1,3); 0, 60/255, 0], 1);
% Getting parameter values
rel_win = params.relative_window; % [-1, 1]*0.8;
del_win = params.delay_window;
bin_size = params.bin_size; % 10*m;
Nb = params.Ns;
% Computing time axis for trials
trial_tx = (rel_win(1) + bin_size/2):bin_size:(rel_win(2) - bin_size/2);

brWin_aux = brWin;
if brWin(1) < 0.12
    brWin_aux = brWin + 0.1;
end
brWin = [brWin_aux; repmat( brWin, Nb-1, 1 )];


%% Organising figures in subfolders
vwKey = sprintf("V%.2f - %.2f s", rel_win);
dwKey = sprintf("D%.2f - %.2f ms", del_win*k);
rwKey = sprintf("R%.2f - %.2f ms", brWin(end,:)*k);
bsKey = sprintf("B%.2f ms", bin_size*k);
subFig = "Beh %s %s";
% Configuration subfolder
% subfigDir = fullfile(figureDir, sprintf(subFig, vwKey, rwKey));
% metaNameFlag = false;
% if exist(subfigDir, "dir")
%     figureDir = subfigDir;
% else
%     if ~mkdir(subfigDir)
%         % Print metadata on figure name
%         metaNameFlag = true;
%         if verbose
%             fprintf(1, "Error while creating subfolder!\n")
%             fprintf(1, "Placing the figures in 'Figure' directory.\n")
%         end
%     else
%         figureDir = subfigDir;
%     end
% end
vec2tr = @(x) reshape( x, params.Nb, [], params.Ns );

% Reconstructing behaviour (Laser OFF only)
mdl_mu = squeeze( mean( mdlAll_ind, 2 ) );
rb = cellfun(@(x) x * mdl_mu, DX([2,3]), fnOpts{:} );
% Reorganising to have trials
rbStack = cellfun(vec2tr, rb, fnOpts{:} );
obStack = cellfun(vec2tr, DX([1,4]), fnOpts{:} );
[~, Nr, Nb] = size( rbStack{1} );
stk = cat( 1, rbStack(1), obStack(1) );
mvpt = zeros( Nb, Nr, 2 );
up_pc = 1.15;
for cs = 1:2
    % Getting maximum speed per trial
    mvpt_aux = arrayfun(@(b) getMaxAbsPerTrial( squeeze( stk{cs}(:,:,b) ), ...
        brWin(b,:), bsWin, trial_tx ), 1:Nb, fnOpts{:} );
    mvpt(:,:,cs) = cat( 1, mvpt_aux{:} );
end
% Normalising each
mvps = max( mvpt, [], 2 ) * up_pc;
ai = squeeze( mean( mvpt ./ mvps, 2 ) );

[f, z_axis, poly_coords] = createPolarPlotPolygons( ai );

[pchObj, ax] = plotPolygons( poly_coords, f, 'clrMap', clrMap );

ai_tot = arrayfun(@(x) area( polyshape( x.Vertices ) ), pchObj );
polLegend = {'Reconstructed', 'Observed'};
polLegend = arrayfun(@(x,y) sprintf([x{:},' %.3f'], y), ...
    polLegend(:), ai_tot(:), fnOpts{:} );

legend( ax, pchObj, polLegend, lgOpts{:} );

arrayfun(@(v,b,y) text( ax, real( z_axis(v) ), ...
    imag( z_axis(v) ), b, vaxOpts{:}, y ), 1:Nb, bodypart_names, ...
    (180*angle( transp( z_axis ) )/pi) - 90 )

results = struct('AmplitudeIndex_pbp', ai, 'AmplitudeIndex', ai_tot, ...
    'AI_perCond', string(polLegend) );
set( f, 'UserData', results )

end