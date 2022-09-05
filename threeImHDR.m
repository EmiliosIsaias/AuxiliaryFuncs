weight_dist = makedist('Normal', 'mu', 2^15, 'sigma', 2^13.5);
weight_fun = @(x) pdf(weight_dist, x);
hdr_img = zeros(size(middle));
wim = hdr_img;
for y = 1:size(middle,1)
    for x = 1:size(middle,2)
        wv = [weight_fun(norm(double(squeeze(bright(y,x,:)))));...
            weight_fun(norm(double(squeeze(middle(y,x,:))))); ...
            weight_fun(norm(double(squeeze(dark(y,x,:)))))]; 
        wv = wv/sum(wv);
        wim = wv;
%         hdr_img(y,x,:) = uint16(double([reshape(bright(y,x,:), 1, 3); ...
%             reshape(middle(y,x,:), 1, 3); ...
%             reshape(dark(y,x,:), 1, 3)])*wv);
    end
end
figure; subplot(1,3,1); imagesc(wim(:,:,1))
subplot(1,3,2); imagesc(wim(:,:,2))
subplot(1,3,3); imagesc(wim(:,:,3))
