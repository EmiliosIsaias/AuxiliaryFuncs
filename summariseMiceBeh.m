function [behTable, miceStruct] = summariseMiceBeh(miceStruct)

fnOpts = {'UniformOutput', false};
dlyPttrn = 'delay\ \d?.\d+\ s$';
dlfPttrn = [dlyPttrn(1:end-1), '\ \+\ L\d+.\d?'];
rxOpts = {'once', 'match', 'freespacing'};

snglFlag = arrayfun(@(m) arrayfun(@(s) string(s.Type) == "single", ...
    m.Sessions), miceStruct, fnOpts{:});

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

[expTypes, ~, expTypeMbrshp] = unique([miceStruct.Structure]);
snglFlag = arrayfun(@(m) arrayfun(@(s) string(s.Type) == "single", ...
    m.Sessions), miceStruct, fnOpts{:});

allCondNames = cellfun(@(m) cellfun(@(t) {t.Conditions}, m), ...
    behTable, fnOpts{:});

allCondNames_aux = cellfun(@(x) cat(1, x{:}), allCondNames, fnOpts{:});
allCondNames_aux = cat(1, allCondNames_aux{:});
[uallCondNames_aux, ~, acnSubs] = unique(allCondNames_aux);
dels = arrayfun(@(x) textscan(x, 'Delay %.3f s'), uallCondNames_aux);
[udels, ~, udSubs] = uniquetol([dels{:}], 0.1 / max([dels{:}]));


behTable = arrayfun(@(m) arrayfun(@(s) {s.DataTable}, m.Sessions), ...
    miceStruct, fnOpts{:});

end