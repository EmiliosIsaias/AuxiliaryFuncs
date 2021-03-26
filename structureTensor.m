function [vMat, lMat, T, nf] = structureTensor(img)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
if ~strcmpi(class(img),'double')
    img = double(img);
end
[Nr, Nc] = size(img);
dx = diff(img,1,2); dx = dx(1:Nr-1,:);
dy = diff(img); dy = dy(:, 1:Nc-1);
nf = cat(3, dx, dy); T = zeros(Nr-1, Nc-1, 4);
for x = 1:2
    for y = 1:2
        T(:,:,(x-1)*2 + y) = nf(:,:,x).*nf(:,:,y);
    end
end
lMat = zeros(Nr-1, Nc-1, 2);

vMat = zeros(Nr-1, Nc-1, 4);
for cc = 1:Nc-1
    for cr = 1:Nr-1
        ct = reshape(T(cr,cc,:),2,2);
        [eV, eAux] = eig(ct);
        lMat(cr,cc,:) = reshape(diag(eAux), 1, 1, 2);
        vMat(cr,cc,:) = reshape(eV, 1, 1, 4);
    end
end
end