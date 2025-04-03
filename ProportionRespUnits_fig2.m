
expandName = @(x) fullfile( x.folder, x.name );
fnOpts = {'UniformOutput', false};
cellcat = @(x,d) cat( d, x{:} );
getMI = @(x,d) diff(x, 1, d) ./ ...
    sum(x, d).*(sum(x,d)>0) + (1.*(sum(x,d)==0 | sum(x,d)< 1e-12));
%%
roller_path = "Z:\Emilio\SuperiorColliculusExperiments\Roller";
% Responsive units in all experiments between 20 and 200 ms
res_paths = dir( fullfile( roller_path, "Batch*", "*", "*", "*", "ephys*", ...
    "Results", "Res VW* ms RW20.00-200.00 ms* PuffAll.mat" ) );
% map_paths = dir( fullfile( roller_path, "Batch*", "*", "*", "*", ...
%     "ephys*", "Results", "Map *.mat" ) );
Nexp = size( res_paths, 1 );

clInfo = cellfun(@(x) getClusterInfo( expandName( ...
    dir( fullfile( getParentDir( x, 1 ), 'cluster_info.tsv' ) ) ) ), ...
    {res_paths.folder}', fnOpts{:} );

rsVars2load = {'Results', 'gclID', 'Counts'};
mfVars2load = {'keyCell', 'resMap'};
%%
% map_res_flag = ismember({map_paths.folder}', {res_paths.folder}');
% mr_subs = find( map_res_flag );

[~, unq_subs] = unique( string( {res_paths(:).folder} ) );
rp_flags = unq_subs == (1:Nexp);
rep_pos = ~sum( rp_flags );
rp_flags( :, rep_pos ) = rp_flags( :, find( rep_pos ) - 1 );
Nu = cellfun(@(x) sum( x.ActiveUnit ), clInfo );
% Nuinit = cumsum( [1; Nu(1:end-1) ] );
% Nuend = cumsum( Nu );
% PSTHall_mu = zeros( 700, sum( Nu ) );
% brAll = [];
% PSTHall = cell( Nexp, 1 );
uSig = cell( Nexp, 1 ); uID = uSig;
uMod = uSig;
uMI = uSig;
%%
for ce = rp_flags'
    if any(find(ce)==102)
        fprintf(1,"What's wrong with this session?")
    end
    nms = {res_paths(ce).name}';
    % View, response, spontaneous
    params = cellfun(@(x) str2double( x ), ...
        cellcat( cellfun(@(y) regexp( y, '\d+\.\d+', 'match' ), ...
        nms, fnOpts{:} ), 1 ) );
    params = params .* -[1,1,-1,-1,1,1];
    min_subs = [];
    sel = 1;
    if nnz( ce ) > 1
        params_diff = diff( params, 1, 1 );
        if ~nnz( params_diff(3:end) )
            % Equivalent results
            ce = find( ce, 1, "first" );
            params = params(1,:);
            nms = nms(1);
        else
            % Spontaneous window messed up. Selecting the furthest and
            % equal to response window's length.
            ce = find( ce );
            sWins = diff( params(:,5:6), 1, 2 );
            rWins = diff( params(:,3:4), 1, 2 );
            eq_flag = sWins == max( rWins );
            if all(eq_flag)
                [~, min_subs] = min( params(:,5:6), [], 1 );
                ce = ce(min_subs(1)); %#ok<*FXSET>
                params = params(min_subs(1),:);
                nms = nms(min_subs(1));
                sel = min_subs(1);
            else
                ce = ce(eq_flag);
                params = params(eq_flag);
                nms = nms(eq_flag);
                sel = nms(eq_flag);
            end
        end
    end

    data_dir = getParentDir( res_paths(ce).folder, 1 );
    condStruct = load( expandName( dir( fullfile( data_dir, ...
        "*analysis.mat" ) ) ), "Conditions", "fs" );
    Conditions = condStruct.Conditions; fs = condStruct.fs;
    Ncond = numel( Conditions );
    p_cond = contains( string({Conditions.name}'), 'control puff', 'IgnoreCase', true );

    mpPath = expandName( dir( fullfile( res_paths(ce).folder, "Map *.mat" ) ) );
    mfStruct = load( mpPath, mfVars2load{:} );
    resMap = mfStruct.resMap; keyCell = mfStruct.keyCell;

    rsPath = expandName( res_paths(ce) );
    rsStruct = load( rsPath, rsVars2load{:} );
    Results = rsStruct.Results; gclID = rsStruct.gclID;
    disp( keyCell(:,4) )
    if numel(params)==6
    configFlag = contains( keyCell(:,1), 'RW20.00-200.00' ) & ...
        contains( keyCell(:,4), 'Control Puff' ) & ...
        contains( keyCell(:,2), sprintf( 'SW%.2f-%.2f', params(5:6) ) );
    else
        fprintf(1, 'What?!\n')
    end
    if nnz(configFlag) > 1 && isscalar(ce)
        configFlag = find( configFlag );
        configFlag = configFlag(sel);
    elseif nnz(configFlag) > 1
        fprintf(1, 'Stop here!\n')
    end
    soe_sub = strfind( keyCell(configFlag,4), 'Control Puff' );
    % uSig{ce} = resMap( keyCell{configFlag,:} );
    % uMod{ce} = zeros( Nu(ce), 1 ); uMI{ce} = uMod{ce};
    % uID{ce} = clInfo{ce}{clInfo{ce}.ActiveUnit==1, "cluster_id" };


    if isempty(soe_sub) 
        if Ncond >= 5
            if ( find( p_cond ) / Ncond) < 0.5 || find( p_cond ) == 3
                Counts = rsStruct.Counts(1,:);
                configFlag = cellfun(@(x) ~isempty(x), regexp( {Results.Combination}, ...
                    '1\s1\ssignrank', 'ignorecase' ) );
            else
                Counts = rsStruct.Counts(end,:);
                configFlag = numel(Results);
            end

        end
    else
        if soe_sub{1} == 1
            Counts = rsStruct.Counts(1,:);
            configFlag = cellfun(@(x) ~isempty(x), regexp( {Results.Combination}, ...
                '1\s1\ssignrank', 'ignorecase' ) );
        else
            Counts = rsStruct.Counts(end,:);
            configFlag = numel(Results);
        end
    end
    Counts = mean( cat( 3, Counts{:} ), 2 );
    uSig{ce} = Results(configFlag).Activity(1).Pvalues;
    uMod{ce} = Results(configFlag).Activity(1).Direction;
    uMI{ce} = getMI( Counts, 3 );
    uID{ce} = gclID;

end