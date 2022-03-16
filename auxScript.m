bstFt = zeros(Nti(1),1,"single");
selSbs = zeros(mxSub, Nti(1), "uint8");
figure; axs = axes("NextPlot","add",Color="none");

for ci = 1:Nti(cc)
    iS = setdiff(1:Nti(1), ci);
    M = [itTimes{1}(iS,1), atTimes{1}];
    scatter(axs, itTimes{1}(iS,1), atTimes{1})
    mdl = fit_poly(itTimes{1}(iS,1), atTimes{1},1);
    [n, d] = getHesseLineForm(mdl);
    line(itTimes{1}(:,1), itTimes{1}(:,1).^[1,0]*mdl, "Color", "k")
    yErr = (M*n - d).^2;
    bstFt(ci) = log(mean(yErr));
    selSbs(:,ci) = iS;
end