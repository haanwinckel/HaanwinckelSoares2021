
%Load moments and create estimator object
tfit = [];
tcf = [];
tpol = [];
load('../CleanData/moments_regions.mat');
id_valid = 1;
for id = 1:100
    id
    filename = ['./estimationOutput/cfResults3_' num2str(id) '.mat'];
    if isempty(dir(filename))
        continue;
    end
    load(filename);
    [ table_fit, table_cf, table_pol ] = makeTablesFromModelStats(st);
    tfit(:,:,id_valid) = cell2mat(table_fit(2:end,2:end));
    tcf(:,:,id_valid) = cell2mat(table_cf(2:end,2:end));
    tpol(:,:,id_valid) = cell2mat(table_pol(2:end,2:end));
   id_valid = id_valid + 1;
end
table_fit(2:end,2:end) = num2cell(std(tfit,0,3));
table_cf(2:end,2:end) = num2cell(std(tcf,0,3));
table_pol(2:end,2:end) = num2cell(std(tpol,0,3));

%Add SE of national moments
rowTitles = {'Labor share';'Almeida Carneiro elast';'Frac 500/100'};
table_fit = [table_fit; ...
    [rowTitles num2cell(sqrt(diag(V((end-2):end,(end-2):end)))) {'';'';''}]];

%Save to Excel
arrays = {table_fit;table_cf;table_pol};
for i = 1:3
    saveCellArrayToExcel(arrays{i},'../Output/standard_errors.xls',i);
end