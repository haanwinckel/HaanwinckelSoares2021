
labels = {'rob_low_gamma';'main';'rob_high_gamma';...
    'rob_sigma_inf';'rob_unemp_ins'};

for r = 1:5
    load(['cleanOutput/cfMain_' num2str(r) '.mat']);
    if r == 2
        makeFigures = true;
    else
        makeFigures = false;
    end
    [ table_fit, table_cf, table_pol ] = ...
        makeTablesFromModelStats(st, makeFigures);
    filename = ['../Output/Tables_' labels{r} '.xls'];
    %Include national moments in table_fit
    load(['cleanOutput/estMain_' num2str(r) '.mat']);
    rowTitles = {'Labor share';'Almeida Carneiro elast';'Frac 500/100'};
    table_fit = [table_fit; ...
        [rowTitles num2cell(...
        [results.nationalTargets results.nationalMoments])]];
    saveCellArrayToExcel(table_fit,filename,1);
    saveCellArrayToExcel(table_cf,filename,2);
    saveCellArrayToExcel(table_pol,filename,3);
end
close all