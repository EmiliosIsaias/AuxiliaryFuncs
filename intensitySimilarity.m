[Nr, Nc] = size(A);
nr = 6.5;
ns = ceil(nr);
nrect = 2*ns + 1;
intSim = @(cp, np, im) exp(-((im(cp(2),cp(1)) - im(np).^2) ./ var(im(np))));
cimg = double(A);
intSimImg = zeros(size(A));
for cc = 1:Nc
    x = repmat((cc-ns):(cc+ns), nrect, 1);
    for cr = 1:Nr
        y = repmat(((cr-ns):(cr+ns))', 1, nrect);
        % circle equation & inside the image (& not centre?)
        coords = cat(3, x, y); C = cat(3, cc, cr);
        relIdx = sqrt(sum((coords - C).^2,3)) <= nr & all(coords >= 1,3)...
            & all(coords <= cat(3, Nc, Nr), 3) & any(coords ~= C, 3); 
        ncp = sum(relIdx(:)); cpSubs = coords(repmat(relIdx(:),2,1)); 
        cpLinSubs = sub2ind([Nr,Nc],cpSubs((ncp+1):(ncp*2)),cpSubs(1:ncp));
        intSimImg(cr,cc) = mean(intSim([cc,cr], cpLinSubs, cimg),'omitnan');
    end
end
figure; imagesc(intSimImg); axis image; axis off;