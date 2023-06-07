function [outMat] = extractBIfromTbl(inTables,Ncc)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
outMat = nan(sum(inTables{:,"Conditions"}=="Control Puff"), Ncc);
rcnt = 0;
for cr = 1:size(inTables,1)
    if inTables{cr,"Conditions"} == "Control Puff"
        rcnt = rcnt + 1;
        ccnt = 1;
    else
        if contains(inTables{cr,"Conditions"}, '+ L')
            ccnt = 3;
        else
            ccnt = 2;
        end
    end
    outMat(rcnt, ccnt) = inTables{cr,"BehaviourIndices"};
end
end