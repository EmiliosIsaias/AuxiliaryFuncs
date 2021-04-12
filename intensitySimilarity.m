function [simMask] = intensitySimilarity(img, radius)
if ~strcmpi(class(img),'double')
    img = double(img);
end
img = gpuArray(img);
[Nr, Nc] = size(img);
if ~exist('radius', 'var')
    radius = 6.5;
end
ns = ceil(radius);
nrect = 2*ns + 1;
intSim = @(cp, np, im) exp(-((im(cp(2),cp(1)) - im(np).^2)) /...
    (2*std(im(np))^2))/(std(im(np))*(sqrt(2*pi)));
simMask = zeros(size(img),'gpuArray');
for cc = 1:Nc
    x = repmat((cc-ns):(cc+ns), nrect, 1);
    for cr = 1:Nr
        y = repmat(((cr-ns):(cr+ns))', 1, nrect);
        coords = cat(3, x, y); C = cat(3, cc, cr);
        % circle equation & inside the image (& not centre?)
        relIdx = sqrt(sum((coords - C).^2,3)) <= radius & all(coords >= 1,3)...
            & all(coords <= cat(3, Nc, Nr), 3) & any(coords ~= C, 3);
        ncp = sum(relIdx(:)); cpSubs = coords(repmat(relIdx(:),2,1));
        cpLinSubs = sub2ind([Nr,Nc],cpSubs((ncp+1):(ncp*2)),cpSubs(1:ncp));
        simMask(cr,cc) = mean(intSim([cc,cr], cpLinSubs, img),'omitnan');
    end
end
figure; imagesc(simMask); axis image; axis off;

end