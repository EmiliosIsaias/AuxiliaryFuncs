mcirn = nan(sum(mcirnTables{:,"Conditions"}=="Control Puff"), 3);
rcnt = 1; ccnt = 1;
for cr = 1:size(mcirn,1)
    if mcirnTables{cr,"Conditions"} ~= "Control Puff"
        
    else
        mcirn(rcnt, ccnt) = 
    end
end