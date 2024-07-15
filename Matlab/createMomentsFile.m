clear

diary 'estimationOutput/log_createMoments.txt';
pmeF='../CleanData/pme_2003_regions.csv';
natF='../CleanData/national_2003.csv';
for i = 1:199
    pmeFBS{i,1}=['../CleanData/Bootstrap/pme_2003_regions_' num2str(i) '.csv'];
    natFBS{i,1}=['../CleanData/Bootstrap/national_2003_' num2str(i) '.csv'];
end

[moments, V, pmeDataSummary, bsMoments, bsDataSummaries] = ...
    getMomentsAndV(pmeF,natF,pmeFBS, natFBS);

save('../CleanData/moments_regions.mat');

diary OFF

%Export share of skilled workers by group
tab = [{'Age','Education','groups','by schooling';...
    'groups','0-7','8-10','11+'};...
    [{'16-19';'20-24';'25-29';'30-59'} num2cell(...
    reshape(pmeDataSummary{1}.eta_by_educ,4,3))]];
saveCellArrayToExcel(tab,'../Output/shareSkilledByGroup.xls',1);