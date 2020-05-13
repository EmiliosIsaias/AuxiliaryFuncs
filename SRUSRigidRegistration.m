function [r,t] = SRUSRigidRegistration(source,target)
cardinSource = size(source,1);
cardinTarget = size(target,1);
rusQ = [];
numEle = min(cardinSource,cardinTarget);
smin = min(min(source(:,3)),min(target(:,3)));
smax = max(max(source(:,3)),max(target(:,3)));
numSliEle = zeros(2,range(source(:,3)));
SliEle = zeros(2,range(target(:,3)));
if cardinSource < cardinTarget
    for ss = smin:smax
        numSliEle(1,ss-target(1,3)+1) = sum(source(:,3)==ss);
        numSliEle(2,ss-target(1,3)+1) = ss;
        SliEle(1,ss-target(1,3)+1) = sum(target(:,3)==ss);
        SliEle(2,ss-target(1,3)+1) = ss;
    end
    idxs = round(SliEle(1,:)*(cardinSource/cardinTarget));
    if sum(idxs) > numEle
        exceedEle = sum(idxs) - numEle;
        [~,bigSlices] = sort(idxs,'descend');
        idxs(bigSlices(1)) = idxs(bigSlices(1)) - exceedEle;
    elseif sum(idxs) < numEle
        missingEle = numEle - sum(idxs);
        [~,bigSlices] = sort(idxs,'descend');
        idxs(bigSlices(1)) = idxs(bigSlices(1)) + missingEle;
    end
    for ss = target(1,3):target(end,3)
        aux = target(target(:,3)==ss,:);
        idx = randperm(SliEle(1,ss-target(1,3)+1));
        idx = idx(1:idxs(ss-target(1,3)+1));
        rusQ = [rusQ;aux(idx,:)]; %#ok<AGROW>
    end
end
A = [source(:,1) -source(:,2) ones(numEle,1) zeros(numEle,1);...
    source(:,2) source(:,1) zeros(numEle,1) ones(numEle,1)];
b = [rusQ(:,1);rusQ(:,2)];
R = pinv(A)*b;
r = R(1:2)/norm(R(1:2));
r = r(1) + 1j*r(2);
t = R(3:4);
t = t(1) + 1j*t(2);
end 