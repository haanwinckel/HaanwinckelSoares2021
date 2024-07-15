function [moments, V, pmeDataSummary, bsMoments, bsDataSummaries] = ...
    getMomentsAndV(pmeFile,nationalFile,...
    pmeBootstrapFileList, nationalBootstrapFileList)

    [moments, pmeDataSummary] = getMoments(pmeFile,nationalFile);
    B = length(pmeBootstrapFileList);
    bsMoments = nan(length(moments),B);
    bsDataSummaries = cell(B,1);
    for b = 1:B
        disp(repmat('x',3,100));
        disp([repelem('x ',24) num2str(b,'%04.f') repelem(' x',24)]);
        disp(repmat('x',3,100));
        [bsMoments(:,b), bsDataSummaries{b}] = ...
            getMoments(pmeBootstrapFileList{b},...
            nationalBootstrapFileList{b});
    end
    V = cov(bsMoments');
end

function vector = summary2vector(dataSummary)

    R = length(dataSummary);
    vector = [];
    for r = 1:R
        ds = dataSummary{r};
        v = [ds.thetaq;ds.informal;ds.lwh;ds.wp_formal;...
            ds.wp_fs_6_10; ds.wp_fs_11p;...
            ds.fs_6_10_for;ds.fs_11p_for;...
            ds.fs_6_10_inf;ds.fs_11p_inf];
        vector = [vector; v];
    end
end

function [moments, dataSummary] = getMoments(pmeFile,nationalFile)
    dataSummary = getDataSummary(pmeFile);
    vector = summary2vector(dataSummary);
    T = readtable(nationalFile);
    nationalMoments = [T.laborShare;T.infElast;...
        T.share_500_100];
    moments = [vector;nationalMoments];
end