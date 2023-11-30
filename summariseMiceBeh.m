function behTable = summariseMiceBeh(miceStruct)

fnOpts = {'UniformOutput', false};

Nm = numel(miceStruct);
behTable = arrayfun(@(m) arrayfun(@(s) {s.DataTable}, m.Sessions), ...
    miceStruct, fnOpts{:});
expTypes = unique([miceStruct.Structure]);
snglFlag = arrayfun(@(m) arrayfun(@(s) string(s.Type) == "single", ...
    m.Sessions), miceStruct, fnOpts{:});

for cm = 1:Nm
    Ns = numel(miceStruct(cm).Sessions);
    if all(snglFlag{cm})
        
    else

    end
    for cs = 1:Ns
        dTbl = miceStruct(cm).Sessions(cs).DataTable;
    end
end


behTable = arrayfun(@(m) arrayfun(@(s) {s.DataTable}, m.Sessions), ...
    miceStruct, fnOpts{:});

behTable{3}{2}.Conditions{1}(:)


end