function [xmice, miceStruct] = summariseMiceBeh(miceStruct)

fnOpts = {'UniformOutput', false};
dlyPttrn = 'delay\ \d?.\d+\ s$';
dlfPttrn = [dlyPttrn(1:end-1), '\ \+\ L\d+.\d?'];
rxOpts = {'once', 'freespacing'};
% Removing mice with no sessions: anatomy, failed, other
miceStruct(arrayfun(@(m) numel(m.Sessions) == 0, miceStruct)) = [];

snglFlag = arrayfun(@(m) arrayfun(@(s) string(s.Type) == "single", ...
    m.Sessions), miceStruct, fnOpts{:});

[expTypes, ~, expTypeMbrshp] = unique([miceStruct.Structure]);

multiSess = cellfun(@(x) find(~x), snglFlag, fnOpts{:});
multiMice = find(cellfun(@(x) ~isempty(x), multiSess));

for cm = multiMice(:)' 
    for cs = multiSess{cm}(:)'
        currTbl = miceStruct(cm).Sessions(cs).DataTable;
        sessDate = miceStruct(cm).Sessions(cs).Date;
        sessDepth = miceStruct(cm).Sessions(cs).Depth;
        Nnt = size(currTbl,1); % Number of new tables
        auxTbls = cell(Nnt,1);
        for cc = 1:size(currTbl,1)
            auxVars = cellfun(@(v) v, currTbl{cc,:}, fnOpts{:});
            if size( auxVars{1}, 1 ) ~= size( auxVars{end}, 1 )
                auxVars{end} = auxVars{end}';
            end
            auxTbls{cc} = table( auxVars{:}, ...
                'VariableNames', currTbl.Properties.VariableNames );
        end
        auxC = 1;
        for ces = [cs, (1:Nnt-1)+numel(miceStruct(cm).Sessions)]
            miceStruct(cm).Sessions(ces).Date = sessDate;
            miceStruct(cm).Sessions(ces).DataTable = auxTbls{auxC};
            miceStruct(cm).Sessions(ces).Type = 'single';
            miceStruct(cm).Sessions(ces).Depth = sessDepth;
            auxC = auxC + 1;
        end
        %currCondNames = [currTbl.Conditions{:}];
        %currBI = [currTbl.BehaviourIndices{:}];
        %[uqCondNames, ~, ucnSubs] = unique(currCondNames);
        
    end
end
%% Grouping condition names across all mice
Nm = numel(miceStruct); Nspm = arrayfun(@(m) numel(m.Sessions), miceStruct);
Ncpspm = arrayfun(@(m) arrayfun(@(s) size(s.DataTable,1), ...
    m.Sessions), miceStruct, fnOpts{:});
Nem = sum(expTypeMbrshp(:) == reshape(unique(expTypeMbrshp),1,[]));
% Going by delay value... It won't work for other conditions, I think
behTable = arrayfun(@(m) arrayfun(@(s) {s.DataTable}, m.Sessions), ...
    miceStruct, fnOpts{:});
allCondNames = cellfun(@(m) cellfun(@(t) {t.Conditions}, m), ...
    behTable, fnOpts{:});
asCondNames = cellfun(@(m) cat(1, m{:}), allCondNames, fnOpts{:});
allCondNames_aux = cellfun(@(x) cat(1, x{:}), allCondNames, fnOpts{:});
allCondNames_aux = cat(1, allCondNames_aux{:});
[uallCondNames_aux, ~, acnSubs] = unique( allCondNames_aux(:) );
dels = arrayfun(@(x) textscan(x, 'Delay %.3f s'), uallCondNames_aux);
ctrlpCond = contains(uallCondNames_aux, 'control puff', 'IgnoreCase', true);
ctrllCond = contains(uallCondNames_aux, 'laser', 'IgnoreCase', true);
[udels, ~, udSubs] = uniquetol([dels{:}], 0.02 / max([dels{:}]));
udSubs2 = nan(size(dels)); udSubs2(~cellfun(@isempty, dels)) = udSubs;
udSubs2(ctrllCond) = -1; udSubs2(ctrlpCond) = 0; 
udSubs2(isnan(udSubs2)) = (1:sum( isnan( udSubs2 ) ) )' + ...
    max( udSubs2, [], "omitmissing" );


% Assuming the same frequency for all
frqs = arrayfun(@(x) regexpi(x, {dlfPttrn, 'hz'}, rxOpts{:}), uallCondNames_aux,...
    fnOpts{:}); frqs = cellfun(@(z) any(z, 2), cellfun(@(y) ...
    ~cellfun(@isempty, y), frqs, fnOpts{:})); 
pairedStimSubs = [udSubs2, frqs];
[~, ~, dfSubs] = unique(pairedStimSubs, 'rows');
dfSubs(isnan(udSubs2)) = NaN; % delAndFrq(any(isnan(delAndFrq),2),:) = [];

condMembership = [allCondNames_aux, acnSubs, udSubs2(acnSubs), single(frqs(acnSubs))];

sessNames = cell(Nm, 1);
sessDepth = sessNames;
miceNames = arrayfun(@(x) cat(1, miceStruct(expTypeMbrshp==x).Name), ...
    unique(expTypeMbrshp), fnOpts{:});
miceMatsTP = cell(Nm, 1);
miceMatsAI = miceMatsTP;
miceMatsPT = miceMatsTP;
miceMatsPA = miceMatsPT;
miceConds = miceMatsTP;
idxTh = cumsum(cellfun(@sum, Ncpspm)); cidx = 1;
% Organising the data per conditions taking into account the experiment
% type
for cm = 1:Nm
    mSessNames = string({miceStruct(cm).Sessions.Date})';
    sessNames{cm} = mSessNames;
    mSessDepth = string({miceStruct(cm).Sessions.Depth})';
    sessDepth{cm} = mSessDepth;
    for cs = 1:numel( miceStruct(cm).Sessions )
        if iscell( miceStruct(cm).Sessions(cs).DataTable.Trial_and_Amp_Indices )
            miceStruct(cm).Sessions(cs).DataTable.Trial_and_Amp_Indices = ...
                cat(1, miceStruct(cm).Sessions(cs).DataTable.Trial_and_Amp_Indices{:} );
        end
    end
    % mBI = arrayfun(@(s) s.DataTable.BehaviourIndices, ...
    %     miceStruct(cm).Sessions, fnOpts{:}); mBI = cat(1, mBI{:});
    mTP = arrayfun(@(s) s.DataTable.Trial_and_Amp_Indices(:,1), ...
        miceStruct(cm).Sessions, fnOpts{:}); mTP = cat(1, mTP{:});
    mAI = arrayfun(@(s) s.DataTable.Trial_and_Amp_Indices(:,2), ...
        miceStruct(cm).Sessions, fnOpts{:}); mAI = cat(1, mAI{:});

    % Polygon unfolding for trial proportions
    mPT = arrayfun(@(s) cellfun(@(c) c(:,1)', s.DataTable.PolygonUnfold, ...
        fnOpts{:}), miceStruct(cm).Sessions, fnOpts{:} );
    mPT = cellfun(@(c) cat( 1, c{:} ), mPT, fnOpts{:} ); 
    mPT = cat( 1, mPT{:} );
    % Polygon unfolding for amplitude index
    mPA = arrayfun(@(s) cellfun(@(c) c(:,2)', s.DataTable.PolygonUnfold, ...
        fnOpts{:}), miceStruct(cm).Sessions, fnOpts{:} );
    mPA = cellfun(@(c) cat( 1, c{:} ), mPA, fnOpts{:} );
    mPA = cat( 1, mPA{:} );

    csCondNames = asCondNames{cm};
    mCondClass = dfSubs(acnSubs(cidx:idxTh(cm)));
    cidx = idxTh(cm) + 1;
    [uAux, ~, muc] = unique(mCondClass); Ncpm_as = numel(uAux);
    mouseMat_TP = nan(Ncpm_as, Nspm(cm));
    mouseMat_AI = mouseMat_TP;
    mouseMat_PP = nan( Ncpm_as, Nspm(cm), 8 );
    mouseMat_PA = mouseMat_PP;
    csCondNames2 = strings(Nspm(cm), Ncpm_as); cidx2 = 1;
    for cc = 1:Nspm(cm)
        cidx3 = cidx2:cidx2-1+Ncpspm{cm}(cc);
        sIdx = muc(cidx3);
        mouseMat_TP(sIdx, cc) = mTP(cidx3);
        mouseMat_AI(sIdx, cc) = mAI(cidx3);
        mouseMat_PP(sIdx, cc, :) = mPT(cidx3, :);
        mouseMat_PA(sIdx, cc, :) = mPA(cidx3, :);
        csCondNames2(cc, sIdx) = csCondNames(cidx3);
        cidx2 = cidx2 + Ncpspm{cm}(cc);
    end
    miceMatsTP{cm} = mouseMat_TP';
    miceMatsAI{cm} = mouseMat_AI';
    miceMatsPT{cm} = mouseMat_PP;
    miceMatsPA{cm} = mouseMat_PA;
    miceConds{cm} = csCondNames2;
end

Nc_s = cellfun(@size, miceMatsTP, fnOpts{:}); Nc_s = cat(1, Nc_s{:});
mxC = arrayfun(@(x) max(Nc_s(expTypeMbrshp == x,:),[],1), ...
    unique(expTypeMbrshp), fnOpts{:});

sessNames = arrayfun(@(x) {sessNames(expTypeMbrshp == x)}, ...
    unique(expTypeMbrshp));
sessDepth = arrayfun(@(x) {sessDepth(expTypeMbrshp == x)}, ...
    unique(expTypeMbrshp));

homMiceMat = arrayfun(@(x) ...
    cellfun(@(m) ...
    padarray(m, mxC{x} - size(m), NaN, "post"), ...
    miceMatsTP(expTypeMbrshp == x), fnOpts{:}), ...
    unique(expTypeMbrshp), fnOpts{:});
homMiceMat = cellfun(@(x) cat(3, x{:}), homMiceMat, fnOpts{:});

% Padding the columns and then the rows of the condition names
homMiceConds = arrayfun(@(m,x) ...
    [miceConds{m}, strings([Nc_s(m,1), mxC{x}(2)-Nc_s(m,2)])], ...
    (1:Nm)', expTypeMbrshp, fnOpts{:});
homMiceConds = arrayfun(@(m,x) ...
    [homMiceConds{m}; strings([mxC{x}(1)-Nc_s(m,1), mxC{x}(2)] ) ], ...
    (1:Nm)', expTypeMbrshp(:), fnOpts{:});
homMiceConds = arrayfun(@(x) cat(3, homMiceConds{expTypeMbrshp == x }), ...
    unique(expTypeMbrshp), fnOpts{:});


%% Experimental group
xmice = struct('ExperimentalGroup', cellstr(expTypes(:)), ...
    'DataTable', homMiceMat(:), ...
    'ConditionNames', homMiceConds(:), ...
    'MiceNames', miceNames(:), ...
    'SessionDate', sessNames(:), ...
    'SessionDepth', sessDepth(:) );

end