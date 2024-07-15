function runBootstrapCounterfactuals(id)

    diaryFileName = ['./estimationOutput/cfLog3_' num2str(id) '.txt'];
    delete(diaryFileName);


    load('../CleanData/moments_regions.mat');
    ME = ModelEstimator(bsMoments(:,id),V,bsDataSummaries{id});
    %%% Recover best result from first step
    filename = ['./estimationOutput/bsEst_' num2str(id) '.txt'];
    fhandle = fopen(filename,'rt');
    thisline = fgetl(fhandle);
    gamma = 0.4366;
    best_loss = inf;
    while true
        thisline = fgetl(fhandle);
        if ~ischar(thisline); break; end  %end of file

        lineStartCheck = 'x = ['; 
        beginning = thisline(1:min(...
            length(thisline),length(lineStartCheck)));
        if strcmp(beginning,lineStartCheck)
            last_x = str2num(thisline(6:(end-1)))';
        end
        lineStartCheck = 'Loss = '; 
        beginning = thisline(1:min(...
            length(thisline),length(lineStartCheck)));
        if strcmp(beginning,lineStartCheck)
            lastChar = min(length(lineStartCheck)+6,length(thisline));
            val = str2double(thisline(...
                (length(lineStartCheck)+1):lastChar));
            if val <= best_loss
                best_loss = val;
                best_x = last_x;
            end
        end
    end
    fclose(fhandle);

    ME.setAllParameters(best_x);

diary(diaryFileName);
disp(gamma);
disp(best_x);
paramsToCleanNonzero={'lambda_for';'lambda_inf'};

for i = 1:ME.numRegions
    for m = 1:length(paramsToCleanNonzero)
        ME.regionModels{i}.(paramsToCleanNonzero{m}) = max(ME.regionModels{i}.(paramsToCleanNonzero{m}),0.001);
    end
    ME.regionModels{i}.verbosity = 5;
    ME.regionModels{i}.gamma = gamma;
end

for i = 1:ME.numRegions
    results.eqs{i,1} = ME.regionModels{i}.solveForEquilibrium();
end

fprintf('\n\n\nNOW PROCEEDING TO runCounterfactuals\n\n\n');


rng(1234);
    bootstrapId = id;

dataFile2012 = ['../CleanData/Bootstrap/pme_2012_regions_' ...
    num2str(bootstrapId) '.csv'];

mwChange = 1.612;
rhoChange = log(1.339);
tau2012 = ME.regionModels{1}.tau-0.0063;
bD2012 = ME.regionModels{1}.bD;
bF2012 = ME.regionModels{1}.bF;
bV2012 = [0.268;0.314];

CM = CounterfactualMaker(ME,results,...
    dataFile2012,mwChange,rhoChange,tau2012,bD2012,...
    bF2012,bV2012,bootstrapId);
st = CM.statTable;
pt = CM.paramTable;
save(['estimationOutput/cfResults3_' num2str(id) '.mat'],'st','pt');

diary OFF
end