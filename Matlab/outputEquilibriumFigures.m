load('cleanOutput/estMain_2.mat');

vertAxisMinMax = [-0.2 1.6;-0.2 1.6;-0.2 6; 0 1];

for r = 1:ME.numRegions
    f = ME.regionModels{r}.plotEquilibrium(...
        results.eqs{r},1:5,vertAxisMinMax);
    saveFigure(f,['../Output/fig_eq_' ME.regionNames{r} '.pdf'],5.2,4.2);
end
close all