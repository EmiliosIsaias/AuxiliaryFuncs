
expandName = @(x) fullfile( x.folder, x.name );
fnOpts = {'UniformOutput', false};
roller_path = "Z:\Emilio\SuperiorColliculusExperiments\Roller";
% Responsive units in all experiments between 20 and 200 ms
res_paths = dir( fullfile( roller_path, "Batch*", "*", "*", "*", "ephys*", ...
    "Results", "Res VW* ms RW20.00-200.00 ms* PuffAll.mat" ) );
map_paths = dir( fullfile( roller_path, "Batch*", "*", "*", "*", ...
    "ephys*", "Results", "Map *.mat" ) );
Nexp = size( res_paths, 1 );

clInfo = cellfun(@(x) getClusterInfo( expandName( ...
    dir( fullfile( getParentDir( x, 1 ), 'cluster_info.tsv' ) ) ) ), ...
    {res_paths.folder}', fnOpts{:} );
%%
map_res_flag = ismember({map_paths.folder}', {res_paths.folder}');
Nu = cellfun(@(x) sum( x.ActiveUnit ), clInfo );
% Nuinit = cumsum( [1; Nu(1:end-1) ] );
% Nuend = cumsum( Nu );
% PSTHall_mu = zeros( 700, sum( Nu ) );
% brAll = [];
% PSTHall = cell( Nexp, 1 );
uSig = cell( Nexp, 1 ); uID = uSig;
uMod = uSig;
uMI = uSig;

for ce = 1:Nexp
    data_dir = getParentDir( res_paths(ce).folder, 1 );
    % condStruct = load( expandName( dir( fullfile( data_dir, "*\*analysis.mat" ) ) ), afVars2load{:} );
    % rstPath = dir( fullfile( data_dir, "*", "*RW20.00-200.00*RelSpkTms.mat" ) );
    % brPath = dir( fullfile( data_dir, "*","BehaviourResult*.mat" ) );
    % mfPath = dir( fullfile( data_dir, "ephys*", "Results", "Res VW* RW20.00-200.00 ms SW*PuffAll.mat") );
    mfPath = res_paths(ce).folder;
    
    % brVars2load = 'behRes';
    % if isempty( brPath )
    %     brPath = dir( fullfile( data_dir, "*\Simple summary.mat" ) );
    %     brVars2load = 'summStruct';
    % end
    % brStruct = load( expandName( brPath ), brVars2load );
    % behRes = brStruct.(brVars2load);
    % ctrlSub = ismember( string( {behRes.ConditionName} ), 'Control Puff' );
    % brAll = cat( 1, brAll, behRes(ctrlSub) );
    % if ~isempty(rstPath) || isscalar( rstPath )
    %     rstCont = load( expandName( rstPath ), rstVars2load{:} );
    % else
    %     fprintf(1, 'Either empty or more than 1 file found!\n');
    %     disp( {rstPath.name} )
    %     continue
    % end
    mftype = 1;
    mfVars2load = {'Results', 'gclID', 'Counts'};
    if isempty( mfPath )
        fprintf( 1, 'Res file not found!\n')
        mfPath = dir( fullfile( data_dir, "ephys*", "Results", "Map*.mat") );
        mfVars2load = {'keyCell', 'resMap'};
        if isempty( mfPath )
            fprintf( 1, 'Unable to load unit response information!\n')
            disp( mfPath )
            continue
        end
        mftype = 2;
    end
    mfStruct = load( expandName( mfPath ), mfVars2load{:} );
    if mftype == 1
        Results = mfStruct.Results; gclID = mfStruct.gclID;
        Counts = mfStruct.Counts(1,:); Counts = mean( cat( 3, Counts{:} ), 2 );
        configFlag = cellfun(@(x) ~isempty(x), regexp( {Results.Combination}, ...
            '1\s1\ssignrank', 'ignorecase' ) );
        uSig{ce} = Results(configFlag).Activity(1).Pvalues;
        uMod{ce} = Results(configFlag).Activity(1).Direction;
        uMI{ce} = getMI( Counts, 3 );
        uID{ce} = gclID;
    else
        resMap = mfStruct.resMap; keyCell = mfStruct.keyCell;
        configFlag = contains( keyCell(:,1), 'RW20.00-200.00' ) & ...
            contains( keyCell(:,4), 'Control Puff' );
        if sum( configFlag ) ~= 1
            fprintf( 1, 'Cannot process this keycell... \n')
            disp( keyCell )
        end
        uSig{ce} = resMap( keyCell{configFlag,:} );
        uMod{ce} = zeros( Nu(ce), 1 ); uMI{ce} = uMod{ce};
        uID{ce} = clInfo{ce}{clInfo{ce}.ActiveUnit==1, "cluster_id" };
    end
    % rstStruct = rstCont.relativeSpkTmsStruct;
    % confStruct = rstCont.configStructure;
    % Conditions = condStruct.Conditions;
    % fs = condStruct.fs;
    % if any( confStruct.Viewing_window_s ~= vWin )
    %     confStruct.Viewing_window_s = vWin;
    % end
    % [PSTH, trial_tx, Na] = getPSTH_perU_perT( rstStruct, confStruct );
    % a = Nuinit(ce); b = Nuend(ce);
    % idx = a:b;
    % PSTHall_mu(:,idx) = squeeze( mean( PSTH{1}, 1 ) );
    % PSTHall(ce) = PSTH(1);
end
clearvars PSTH brStruct ctrlSub rstStruct confStruct brPath brVars2load ...
    rstCont a b idx;
ai_pt = arrayfun(@(s) getAIperTrial( s ), brAll, fnOpts{:} );
