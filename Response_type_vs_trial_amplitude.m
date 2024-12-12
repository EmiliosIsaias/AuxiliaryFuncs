expandName = @(x) fullfile( x.folder, x.name );
roller_path = "Z:\Emilio\SuperiorColliculusExperiments\Roller";
owfFlag = false;
m = 1e-3; k = 1e3;
fnOpts = {'UniformOutput', false};
%%
distPercent = 0.25;
pePaths = ["Batch12_ephys.e\MC\vGlut1\221206_C+F_2100",...
    "Batch11_ephys.MC\ChR2\Rb69\221128_F+C_2370", ...
    "Batch17_ephys.MC\eOPN3\Rb27\231103_C_2500", ...
    "Batch10_ephys\MC\WTg64\221111_C+PTX10microM_1900", ...
    "Batch10_ephys\MC\WTg64\221110_C_1900"];
vgePaths = ["Batch18_ephys\MC\GADi41\240213_F+C_2450", ...
    "Batch16_ephys\MC\GADi21\231023_C+F_2350", ...
    "Batch12_ephys.e\MC\vGlut5\221210_F+C_2450", ...
    "Batch11_ephys.MC\eOPN3\WT67\221202_C_2240", ...
    "Batch11_ephys.MC\ChR2\Rb71\221130_F+C_2100", ...
    "Batch11_ephys.MC\ChR2\Rb69\221129_C+F+Mus_2000", ...
    "Batch10_ephys\MC\WTg65\221115_E+PTX20microM_1900", ...
    "Batch8_ephys\MC\GADe55\220825_C+F_DV1750", ...
    "Batch7_ephys\MC\GADi52\220809_C+F_1900"];
aePaths = cat( 1, pePaths', vgePaths' );
rstVars2load = {'relativeSpkTmsStruct', 'configStructure'};
afVars2load = {'Conditions', 'fs'};
%%
vWin = [-300,400]*m;
clInfo = arrayfun(@(x) getClusterInfo( expandName( ...
    dir( fullfile( roller_path, x, '*\cluster_info.tsv' ) ) ) ), ...
    aePaths, fnOpts{:} );
%%
Nu = cellfun(@(x) sum( x.ActiveUnit ), clInfo );
Nuinit = cumsum( [1; Nu(1:end-1) ] );
Nuend = cumsum( Nu );
PSTHall_mu = zeros( 700, sum( Nu ) );
brAll = [];
PSTHall = cell( numel( aePaths ), 1 );
%%
for ce = 1:numel( aePaths )
    data_dir = fullfile( roller_path, aePaths(ce) );
    % condStruct = load( expandName( dir( fullfile( data_dir, "*\*analysis.mat" ) ) ), afVars2load{:} );
    rstPath = dir( fullfile( data_dir, "*\*RW20.00-200.00*RelSpkTms.mat" ) );
    brPath = dir( fullfile( data_dir, "*\BehaviourResult*.mat" ) );
    brVars2load = 'behRes';
    if isempty( brPath )
        brPath = dir( fullfile( data_dir, "*\Simple summary.mat" ) );
        brVars2load = 'summStruct';
    end
    brStruct = load( expandName( brPath ), brVars2load );
    behRes = brStruct.(brVars2load);
    ctrlSub = ismember( string( {behRes.ConditionName} ), 'Control Puff' );
    brAll = cat( 1, brAll, behRes(ctrlSub) );
    if ~isempty(rstPath) || numel( rstPath ) == 1
        rstCont = load( expandName( rstPath ), rstVars2load{:} );
    else
        fprintf(1, 'Either empty or more than 1 file found!\n');
        disp( {rstPath.name} )
        continue
    end
    rstStruct = rstCont.relativeSpkTmsStruct;
    confStruct = rstCont.configStructure;
    % Conditions = condStruct.Conditions;
    % fs = condStruct.fs;
    if any( confStruct.Viewing_window_s ~= vWin )
        confStruct.Viewing_window_s = vWin;
    end
    [PSTH, trial_tx, Na] = getPSTH_perU_perT( rstStruct, confStruct );
    a = Nuinit(ce); b = Nuend(ce);
    idx = a:b;
    PSTHall_mu(:,idx) = squeeze( mean( PSTH{1}, 1 ) );
    PSTHall(ce) = PSTH(1);
end
clearvars PSTH brStruct ctrlSub rstStruct confStruct brPath brVars2load ...
    rstCont a b idx;
ai_pt = arrayfun(@(s) getAIperTrial( s ), brAll, fnOpts{:} );
%%
[coeff, score, ~, ~, explained] = pca( zscore( ...
    PSTHall_mu(my_xor( trial_tx < [0,350]*m),:) ) );

y = pdist( coeff(:,1:3), 'euclidean' );
z = linkage( y, 'ward' );
rtm = cluster( z, 'cutoff', max(z(:,3))*distPercent, 'criterion', 'distance' );
Ncl = max( unique( rtm ) );
rl_mu = arrayfun(@(x) mean( zscore( PSTHall_mu(:,rtm==x) ), 2 ), ...
    unique(rtm), fnOpts{:} );
rl_mu = cat( 2, rl_mu{:} );
%%
respWins = [20,200; % Whole window
    20,50;          % sensory window
    50,200] * m;    % motor window in milliseconds
sponWins = -flip( respWins, 2 );
getBoolWindows = @(w) arrayfun(@(x) my_xor( trial_tx > w(x,:) ), ...
    1:size( w, 1 ), fnOpts{:} );
cellcat = @(d,x) cat( d, x{:} );

respFlags = cellcat( 2, getBoolWindows( respWins ) );
sponFlags = cellcat( 2, getBoolWindows( sponWins ) );
%%
% Response validation
ruFlags = false( sum( Nu ), 3 );
for cp = 1:numel(PSTHall)
    auxP = false( Nu(cp), 3 );
    for cu = 1:Nu(cp)
        for cr = 1:size( respWins, 1 )
            d1 = squeeze( sum( PSTHall{cp}(:,sponFlags(:,cr),cu), 2 ) );
            d2 = squeeze( sum( PSTHall{cp}(:,respFlags(:,cr),cu), 2 ) );
            [~,auxP(cu,cr)] = poisson_means_etest( d1, d2, 0.05 );
        end
    end
    ruFlags(Nuinit(cp):Nuend(cp),:) = auxP;
end

%% 
