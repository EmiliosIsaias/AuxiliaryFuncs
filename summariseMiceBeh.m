function [behTable, miceStruct] = summariseMiceBeh(miceStruct)

fnOpts = {'UniformOutput', false};
dlyPttrn = 'delay\ \d?.\d+\ s$';
dlfPttrn = [dlyPttrn(1:end-1), '\ \+\ L\d+.\d?'];
rxOpts = {'once', 'freespacing'};

snglFlag = arrayfun(@(m) arrayfun(@(s) string(s.Type) == "single", ...
    m.Sessions), miceStruct, fnOpts{:});

[expTypes, ~, expTypeMbrshp] = unique([miceStruct.Structure]);

multiSess = cellfun(@(x) find(~x), snglFlag, fnOpts{:});
multiMice = find(cellfun(@(x) ~isempty(x), multiSess));

for cm = multiMice(:)' 
    for cs = multiSess{cm}(:)'
        currTbl = miceStruct(cm).Sessions(cs).DataTable;
        sessDate = miceStruct(cm).Sessions(cs).Date;
        Nnt = size(currTbl,1); % Number of new tables
        auxTbls = cell(Nnt,1);
        for cc = 1:size(currTbl,1)
            auxVars = cellfun(@(v) v(:), currTbl{cc,:}, fnOpts{:});
            auxTbls{cc} = table(auxVars{:}, ...
                'VariableNames',currTbl.Properties.VariableNames);
        end
        for ces = cs:cs+Nnt-1
            miceStruct(cm).Sessions(ces).Date = sessDate;
            miceStruct(cm).Sessions(ces).DataTable = auxTbls{ces-cs+1};
            miceStruct(cm).Sessions(ces).Type = 'single';
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

behTable = arrayfun(@(m) arrayfun(@(s) {s.DataTable}, m.Sessions), ...
    miceStruct, fnOpts{:});
allCondNames = cellfun(@(m) cellfun(@(t) {t.Conditions}, m), ...
    behTable, fnOpts{:});
asCondNames = cellfun(@(m) cat(1, m{:}), allCondNames, fnOpts{:});
allCondNames_aux = cellfun(@(x) cat(1, x{:}), allCondNames, fnOpts{:});
allCondNames_aux = cat(1, allCondNames_aux{:});
[uallCondNames_aux, ~, acnSubs] = unique(allCondNames_aux(:));
dels = arrayfun(@(x) textscan(x, 'Delay %.3f s'), uallCondNames_aux);
ctrlCond = contains(uallCondNames_aux, 'control', 'IgnoreCase', true);
[udels, ~, udSubs] = uniquetol([dels{:}], 0.1 / max([dels{:}]));
udSubs2 = nan(size(dels)); udSubs2(~cellfun(@isempty, dels)) = udSubs;
udSubs2(ctrlCond) = 0; udSubs2(isnan(udSubs2)) = find(isnan(udSubs2));


% Assuming the same frequency for all
frqs = arrayfun(@(x) regexpi(x, dlfPttrn, rxOpts{:}), uallCondNames_aux,...
    fnOpts{:}); frqs = ~cellfun(@isempty, frqs); 
frqSubs = find(frqs, 1, "first");
pairedStimSubs = [udSubs2, frqs];
[~, ~, dfSubs] = unique(pairedStimSubs, 'rows');
dfSubs(isnan(udSubs2)) = NaN; % delAndFrq(any(isnan(delAndFrq),2),:) = [];

condMembership = [allCondNames_aux, acnSubs, udSubs2(acnSubs), single(frqs(acnSubs))];

sessNames = cell(Nm, 1);
condNames = cell(numel(expTypes), 1);
miceNames = cat(1, miceStruct.Name);
miceMats = cell(Nm, 1);
idxTh = cumsum(cellfun(@sum, Ncpspm)); cidx = 1;
for cm = 1:Nm
    mSessNames = string({miceStruct(cm).Sessions.Date})';
    sessNames{cm} = mSessNames;
    mBI = arrayfun(@(s) s.DataTable.BehaviourIndices, ...
        miceStruct(cm).Sessions, fnOpts{:}); mBI = cat(1, mBI{:});
    mCondClass = dfSubs(acnSubs(cidx:idxTh(cm)));
    cidx = idxTh(cm) + 1;
    [uAux, ~, muc] = unique(mCondClass); Ncpm_as = numel(uAux);
    mouseMat = nan(Ncpm_as, Nspm(cm));
    cidx2 = 1; 
    for cc = 1:Nspm(cm)
        cidx3 = cidx2:cidx2-1+Ncpspm{cm}(cc);
        sIdx = muc(cidx3);
        mouseMat(sIdx, cc) = mBI(cidx3);
        cidx2 = cidx2 + Ncpspm{cm}(cc);
    end
    miceMats{cm} = mouseMat';
end

Nc_s = cellfun(@size, miceMats, fnOpts{:});
mxC = arrayfun(@(x) max(cat(1, Nc_s{expTypeMbrshp == x}),[],1), ...
    unique(expTypeMbrshp), fnOpts{:});

%% Experimental group
xmice = struct('ExperimentalGroup', {expTypes{:}}, ...
    'DataTable', [], ...
    'MiceNames', arrayfun(@(x) cat(1, miceNames(expTypeMbrshp == x)), ...
    1:numel(expTypes), fnOpts{:}), ...
    'SessionDate', []);


for ce = 1:numel(expTypes)
    xmice(ce).DataTable = nan(sum(expTypeMbrshp == ce), mxC{ce}(1), mxC{ce}(2));

end

end