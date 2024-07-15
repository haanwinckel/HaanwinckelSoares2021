function makeCounterfactuals(resultsFile,outputFile, logFile, bootstrapId, runRobustnessNoEduc)

rng(1234);

if nargin < 4
    bootstrapId = 0;
end

diary(logFile);

load(resultsFile);

if ME.numRegions == 1
    dataFile2012 = '../CleanData/pme_2012_national.csv';
elseif bootstrapId == 0
    dataFile2012 = '../CleanData/pme_2012_regions.csv';
else
    dataFile2012 = ['../CleanData/Bootstrap/pme_2012_regions_' ...
        num2str(bootstrapId) '.csv'];
end
mwChange = 1.612;
rhoChange = log(1.339);
tau2012 = ME.regionModels{1}.tau-0.0063;
bD2012 = ME.regionModels{1}.bD;
bF2012 = ME.regionModels{1}.bF;
bV2012 = [0.268;0.314];

CM = CounterfactualMaker(ME,results,...
    dataFile2012,mwChange,rhoChange,tau2012,bD2012,...
    bF2012,bV2012,bootstrapId, runRobustnessNoEduc);
st = CM.statTable;
pt = CM.paramTable;
save(outputFile,'st','pt');

diary OFF