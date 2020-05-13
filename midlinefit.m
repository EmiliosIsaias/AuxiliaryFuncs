function midline = midlinefit(pts)

err_old = 0;
iter = 1;
count = 0;
cMdl_old = [NaN;NaN];
% Figure
while count < 512 && iter <= 2046
    rIdx = randi(length(pts),round(0.2*length(pts)),1);
    mData = [pts(rIdx,1), pts(rIdx,2)];
    [cMdl,~,err_new]=fit_poly(mData(:,1),mData(:,2),1);
    % err_aux = lineerror(cMdl,pts);
    if err_old < err_new
        % Plot the lines between the estimations
        cMdl_old = cMdl;
        err_old = err_new;
    else
        count = count + 1;
    end
    iter = iter + 1;
end
midline = cMdl_old;
% hold on;plot(pts(:,1),pts(:,1)*cMdl_old(1)+cMdl_old(2),'--r')
end
% figure();plot(pts(:,1),pts(:,2),'LineStyle','none','Marker','o');
 % hold on;plot(pts(:,1),pts(:,1)*cMdl(1)+cMdl(2));text(pts(end,1),...
        %     pts(end,1)*cMdl(1)+cMdl(2),num2str(err_new))