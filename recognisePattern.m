function [mdl] = recognisePattern(signals, timeAx)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
fnOpts = {'UniformOutput', false};
er = 1;
o = 1;
thld = 1/3;
hoFlag = true(size(signals,2),1);
mdl = cell(size(signals,2),1);
ssme = zeros(size(signals,2),1, 'single');
while er > thld && o <= 4
    mdl(hoFlag) = arrayfun(@(t) fit_poly(timeAx(timeAx<0), ...
        signals(timeAx<0,t), o), find(hoFlag), fnOpts{:});
    ssme(hoFlag) = arrayfun(@(t) sqrt(sum((signals(timeAx<0,t) - ...
        (timeAx(timeAx(:)<0).^(o:-1:0) * mdl{t})).^2,1)), find(hoFlag));
    hoFlag(hoFlag) = ssme(hoFlag) > thld;
    o = o + 1;
    er = mean(ssme);
end

end