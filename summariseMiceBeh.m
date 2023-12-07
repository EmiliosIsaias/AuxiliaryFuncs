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

Nm = numel(miceStruct); Nspm = arrayfun(@(m) numel(m.Sessions), miceStruct);
Ncpspm = arrayfun(@(m) arrayfun(@(s) size(s.DataTable,1), ...
    m.Sessions), miceStruct, fnOpts{:});

behTable = arrayfun(@(m) arrayfun(@(s) {s.DataTable}, m.Sessions), ...
    miceStruct, fnOpts{:});
allCondNames = cellfun(@(m) cellfun(@(t) {t.Conditions}, m), ...
    behTable, fnOpts{:});

allCondNames_aux = cellfun(@(x) cat(1, x{:}), allCondNames, fnOpts{:});
allCondNames_aux = cat(1, allCondNames_aux{:});
[uallCondNames_aux, ~, acnSubs] = unique(allCondNames_aux(:));
dels = arrayfun(@(x) textscan(x, 'Delay %.3f s'), uallCondNames_aux);
ctrlCond = contains(uallCondNames_aux, 'control', 'IgnoreCase', true);
[udels, ~, udSubs] = uniquetol([dels{:}], 0.1 / max([dels{:}]));
udSubs2 = nan(size(dels)); udSubs2(~cellfun(@isempty, dels)) = udSubs;
udSubs2(ctrlCond) = 0;

condMembership = [allCondNames_aux, acnSubs, udSubs2(acnSubs)];
% Assuming the same frequency for all
frqs = arrayfun(@(x) regexpi(x, dlfPttrn, rxOpts{:}), uallCondNames_aux,...
    fnOpts{:}); frqs = ~cellfun(@isempty, frqs); 
frqSubs = find(frqs, 1, "first");
pairedStimSubs = [udSubs2, frqs];
[delAndFrq, ~, dfSubs] = unique(pairedStimSubs, 'rows');
dfSubs(isnan(udSubs2)) = NaN; delAndFrq(any(isnan(delAndFrq),2),:) = [];


%% Experimental group
idxTh = cumsum(cellfun(@sum, Ncpspm));
mscSubs = zeros(sum(cellfun(@sum,Ncpspm)), 4, "single");
xmice = struct('ExperimentalGroup', {expTypes{:}}, ...
    'DataTable', [], ...
    'MiceNames',[], ...
    'SessionDate', []);
cidx = 1;
for ce = 1:numel(expTypes)
    % Experimental group
    cx = find(expTypeMbrshp == ce, 1,'last');
    exIdx = acnSubs(cidx:idxTh(cx));
    etpcIdx = pairedStimSubs(exIdx,:);
    condType = dfSubs(exIdx);
    

    cidx = idxTh(cx) + 1;
end
% %%
% idxTh = cumsum(cellfun(@sum, Ncpspm));
% mscSubs = zeros(sum(cellfun(@sum,Ncpspm)), 4, "single");
% xmice = struct('ExperimentalGroup', {expTypes{:}}, ...
%     'DataTable', [], ...
%     'MiceNames',[], ...
%     'SessionDate', []);
% cidx = 1;
% for cm = 1:Nm
%     % Mouse
%     for cs = 1:Nspm(cm)
%         % Session
%         pairedStimSubs(acnSubs(cidx:idxTh(find(expTypeMbrshp==1,1,'last'))),:)
%         udSubs2(acnSubs(cidx:idxTh(cm)))
%         dfSubs(acnSubs(cidx:idxTh(cm)))
%         allCondNames_aux(cidx:idxTh(cm))
%         for cc = 1:Ncpspm{cm}(cs)
%             % Condition
%             
%             mscSubs(idx,:) = [expTypeMbrshp(cm), cm, cs, cc];
%             idx = idx + 1;
%         end
%     end
% end
end