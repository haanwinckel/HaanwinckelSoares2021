
%Load moments and create estimator object
paramTable = [];
load('../CleanData/moments_regions.mat');
for id = 0:100
    id
    if id == 0 %Baseline model
        load('cleanOutput/estMain_2.mat');
    else %Create a ModelEstimator object with bootstrapped parameters
        ME = ModelEstimator(bsMoments(:,id),V,bsDataSummaries{id});
        %%% Recover best result from first step
        filename = ['./estimationOutput/bsEst_' num2str(id) '.txt'];
        if isempty(dir(filename))
            return;
        end
        fhandle = fopen(filename,'rt');
        thisline = fgetl(fhandle);
        gamma = str2double(thisline((end-5):end));
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
    end
    %Store parameters
    paramNames=[...
                ME.commonParamNames;...
                ME.regionSpecificParamNames;...
                {'lambda_for';'lambda_inf'}];
    for i = 1:ME.numRegions
        m = ME.regionModels{i};
        bstrapId = repmat(id,2,1);
        region = {ME.regionNames{i};ME.regionNames{i}};
        skillLevel = {'Skilled';'Unskilled'};
        skillShare = [m.eta;1-m.eta];                
        T = table(bstrapId,region,skillLevel,skillShare);
        for p = 1:length(paramNames)
            pn = paramNames{p};
            val = m.(pn);
            if length(val) == 1
                val = [val;val];
            end
            T.(pn) = val;
        end
        paramTable = [paramTable; T];
    end
end
paramTable.PD=log(paramTable.rhoGridMax2Min)/5;
paramTable.rhoGridMax2Min=[];
paramTable.Pr=paramTable.coefRho1_mu;
paramTable.coefRho1_mu=[];

%National means of parameters
load('cleanOutput/estMain_2.mat');
weightedMean=@(x) sum(x.*ME.regionSizes);
paramTableNational=varfun(weightedMean,paramTable(:,[1 3 4 13:17 19]),'GroupingVariables',{'bstrapId','skillLevel'});

%Make output tables
TO = cell(4,1);

%Estimated params
TO{1} = paramTable(paramTable.bstrapId == 0,:);
%Estimated params, national
TO{2} = paramTableNational(paramTableNational.bstrapId == 0,:);

%Standard errors, all
paramTable = paramTable(paramTable.bstrapId~=0,:);
TO{3}=grpstats(paramTable,{'region','skillLevel'},'std');
TO{3}.std_bstrapId=[];
TO{3}.GroupCount=[];

%Standard errors, national means
paramTableNational = paramTableNational(paramTableNational.bstrapId~=0,:);
TO{4}=grpstats(paramTableNational,{'skillLevel'},'std');
TO{4}.std_bstrapId=[];
TO{4}.GroupCount=[];

%Export
filename = '../Output/parameters.xls';
for i = 1:4
    writetable(TO{i},filename,'sheet',i);
end