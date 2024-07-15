classdef ModelEstimator < handle
    properties (SetAccess = private, GetAccess = public)
        numRegions;
        regionNames;
        regionCodes;
        dataSummary;
        regionModels;
        targetMoments;
        momentsCovMatrix;
        regionSizes;
    end
    properties (SetAccess = public, GetAccess = public)
        commonParamNames = {'alpha';'coefRho2';'B_exp';'B_z_ref';'xi_cons';
            'xi_s';'xi_rel_s';'xi_rel_inf';'rhoGridMax2Min'};
        regionSpecificParamNames = {'A';'T';'coefRho1_mu';'ut_unemp'};         
        numericalOptions_fminsearch = optimset(...
            'display','iter','maxiter',1000,'maxfunevals',10000,...
            'tolx',1e-5,'tolfun',1e-6);        
        numericalOptions_fminunc = optimoptions('fminunc',...
            'display','none','algorithm','quasi-newton',...
            'tolx',1e-5,'tolfun',1e-6,'maxiter',10000,'maxfunevals',10000,...
            'SpecifyObjectiveGradient',true);
        finiteDiffDelta = 1e-5;
        initialSampleDraws = 10;  
        rangeSampleX = 1;
        messageFile = 'screen';     
        centerX = [];
        weightOption = ModelEstimator.sizeWeights;
        algOption = ModelEstimator.alg_fminunc;
    end    
    properties (Constant) %Named options for clarity of code
        optimalWeights = 1;
        sizeWeights = 2;
        alg_fminsearch = 1;
        alg_fminunc = 2;
    end
    
    methods (Access = public)
        
        function this = ModelEstimator(targetMoments,...
                momentsCovMatrix,dataSummary)
            this.numRegions = length(dataSummary);
            multiRegion = this.numRegions > 1;
            if multiRegion
                this.regionNames = {'Recife';'Salvador';'Belo Horizonte';...
                    'Rio de Janeiro';'Sao Paulo';'Porto Alegre'};
                this.regionCodes = [26;29;31;33;35;43];
                for i = 1:this.numRegions
                    this.regionModels{i,1} = Model();
                    this.regionModels{i,1}.outputLabel = this.regionNames{i};
                    this.regionSizes(i,1) = dataSummary{i}.regionSize;
                end
                this.regionSizes = this.regionSizes/sum(this.regionSizes);
            else 
                this.regionNames = {'Brazil'};
                this.regionCodes = [0];
                this.regionModels{1,1} = Model();
                this.regionModels{1,1}.outputLabel = this.regionNames{1};
                this.regionSizes = 1;
            end
            
            this.setTargetsAndObservedRegionalParams(targetMoments,...
                momentsCovMatrix,dataSummary);
        end
        
        function estimateObject = run(this)
            this.initializeMessageFile();
            this.printMessage(['DRAWING INITIAL SAMPLE\nCurrent time:\n'...
                datestr(datetime('now')) '\n'],true);
            
            if this.initialSampleDraws > 1
                sample_x = this.sampleInitialPoints();
            elseif isempty(this.centerX)
                sample_x = this.getAllParameters();
            else
                sample_x = this.centerX;
            end
            numValidSP = size(sample_x,2);
            
            this.printMessage(['STARTING ESTIMATION\nCurrent time:\n'...
                datestr(datetime('now')) '\n'],true);
            
            succeeded = false;
            attempt = 0;
            while ~succeeded && attempt < numValidSP
                attempt = attempt + 1;
                try                    
                    this.printMessage(['\n\n*** SAMPLE #' num2str(attempt) '/' ...
                        num2str(this.initialSampleDraws) ' ***\n']);
                    ini_x = sample_x(:,attempt); 
                    if this.algOption == ModelEstimator.alg_fminsearch
                        [opt_x,loss,flag,optim_output] = fminsearch(...
                            @(x) this.estimationObjFun(x,false),ini_x,...
                            this.numericalOptions_fminsearch);
                    elseif this.algOption == ModelEstimator.alg_fminunc                                         
                        [opt_x,loss,flag,optim_output] = fminunc(...
                            @(x) this.estimationObjFun(x,true),ini_x,...
                            this.numericalOptions_fminunc);
                    else
                        error('Unrecognized algorithm option.');
                    end
                    this.printMessage(['\n\nOptimization complete, flag=' ...
                        num2str(flag) '\nCurrent time:\n'...
                        datestr(datetime('now')) '\nOptimal point:']); 

                    [regionalMoments, nationalContribs, eqs] = ...
                        this.evaluateAllRegions(opt_x);         
                    loss = this.lossVal(regionalMoments,nationalContribs,true,opt_x);  
                    succeeded = true;
                catch excp
                    this.printMessage(['\nAttempt has failed.\nCurrent time:\n'...
                        datestr(datetime('now')) '\n' getReport(excp)],true);                    
                end
            end
            if succeeded  
                this.printMessage(['\nExiting: success.\nCurrent time:\n'...
                    datestr(datetime('now'))],true);               
                estimateObject.x = opt_x;
                estimateObject.eqs = eqs;
                estimateObject.models = this.regionModels;
                estimateObject.flag = flag;
                estimateObject.optim_output = optim_output;
                estimateObject.loss = loss;
            else
                this.printMessage(['\nExiting: failed.\nCurrent time:\n'...
                    datestr(datetime('now'))],true);  
                estimateObject = [];
            end
        end
        
        function setTargetsAndObservedRegionalParams(this,targetMoments,...
                momentsCovMatrix,dataSummary)     
            [this.targetMoments, this.momentsCovMatrix, this.dataSummary]...
                = deal(targetMoments,momentsCovMatrix,dataSummary);            
            for iRegion = 1:this.numRegions
                mom = this.dataSummary{iRegion};
                mod = this.regionModels{iRegion};
                this.dataSummary{iRegion}.q = ...
                    mom.thetaq.^(mod.E/(mod.E-1))...
                    .*(mod.D).^(-1/(mod.E-1));                
                mod.lambda_for = mom.lambda_for;
                mod.lambda_inf = mom.lambda_inf;
                mod.eta = mom.eta;
            end
        end      
        
        function keepOnlySomeRegions(this,whichRegions)
            regionsKept = ismember((1:this.numRegions)',whichRegions);
            this.numRegions=length(whichRegions);
            this.regionNames=this.regionNames(whichRegions);
            this.regionCodes=this.regionCodes(whichRegions);
            this.dataSummary=this.dataSummary(whichRegions);
            this.regionModels=this.regionModels(whichRegions);
            momentsToKeep = [repelem(regionsKept,20,1);true(3,1)];
            this.targetMoments = this.targetMoments(momentsToKeep);
            this.momentsCovMatrix = this.momentsCovMatrix(...
                momentsToKeep,momentsToKeep);
            this.regionSizes=this.regionSizes(whichRegions);
            this.regionSizes=this.regionSizes/sum(this.regionSizes);
        end
    end
    
    methods(Access = public) %Private
        
        function [sample_x, loss] = sampleInitialPoints(this)
            %Make equilibrium computation quicker and less robust, just for
            %drawing the sample
            for r = 1:this.numRegions
                opt_equilibrium{r,1} = this.regionModels{r}.opt_equilibrium;       
                this.regionModels{r}.opt_equilibrium = optimoptions(...
                    this.regionModels{r}.opt_equilibrium,...
                    'tolx',1e-6,'tolfun',1e-6,'display','none',...
                    'SpecifyObjectiveGradient',true,'MaxIter',50);         
                equilibriumMaxAttempts{r,1} = ...
                    this.regionModels{r}.equilibriumMaxAttempts;
                this.regionModels{r}.equilibriumMaxAttempts = 1;
                max_error_equilibrium{r,1} = ...
                    this.regionModels{r}.max_error_equilibrium;
                this.regionModels{r}.max_error_equilibrium = 1e-3;
            end
            if isempty(this.centerX)
                this.centerX = this.getAllParameters();
            end
            random_shifters = this.rangeSampleX * ...
                (rand(length(this.centerX),this.initialSampleDraws)-0.5);
            sample_x = this.centerX + random_shifters;
            loss = nan(this.initialSampleDraws,1);
            for i = 1:this.initialSampleDraws
                this.printMessage(['\n\n*** SAMPLE #' num2str(i) '/' ...
                    num2str(this.initialSampleDraws) ' ***\n']);
                try
                    [regionalMoments, nationalContribs] = ...
                        this.evaluateAllRegions(sample_x(:,i));
                    loss(i) = this.lossVal(regionalMoments,...
                        nationalContribs,true,sample_x(:,i));
                catch excp
                    this.printMessage(['Error calculating loss function.\nx='...
                        num2str(sample_x(:,i)',16) '\n' getReport(excp)]);
                    loss(i) = inf;
                end
            end
            invalid = isnan(loss) | isinf(loss);
            sample_x(:,invalid) = [];
            loss(invalid) = [];
            [loss, order] = sort(loss);
            sample_x = sample_x(:,order);
            %Restore previous optimization properties            
            for r = 1:this.numRegions
                this.regionModels{r}.opt_equilibrium = opt_equilibrium{r,1};        
                this.regionModels{r}.equilibriumMaxAttempts = ...
                    equilibriumMaxAttempts{r,1};
                this.regionModels{r}.max_error_equilibrium = ...
                    max_error_equilibrium{r,1};
            end
        end
        
        function [predictedRegionalMoments,predictedNationalContributions, eq] = lossVector(this,iRegion,referenceEq)
            if nargin < 3 || isempty(referenceEq)
                ini_rU = 0.6*this.dataSummary{iRegion}.wh;
                ini_q = this.dataSummary{iRegion}.q;
            else
                ini_rU = referenceEq.rU;
                ini_q = referenceEq.q;
            end
            eq = this.regionModels{iRegion}.solveForEquilibrium(...
                ini_rU,ini_q,referenceEq);
            if isempty(eq)
                error('ModelEstimator:eqFailed',...
                    'Problem calculating lossVector.');
            end
            ms = this.regionModels{iRegion}.clone();
            delta = 1e-4;
            ms.coefRho1_mu = ms.coefRho1_mu+delta;
            [~,~,eqs] = ms.equilibriumConditions(eq.rU,eq.q);
            if isempty(eqs)
                error('ModelEstimator:eqFailed',...
                    'Problem calculating d inf / d log coefRho1.');
            end
            eqs = this.regionModels{iRegion}.addEquilibriumStatistics(eqs);
            inf_elast = (eqs.informality-eq.informality)/delta;
            predictedRegionalMoments = [...
                eq.theta.*eq.q;...
                eq.informality_by_skill;...
                eq.mean_lw;...
                eq.wp_formal;...
                eq.wp_fs_6_10;...
                eq.wp_fs_11p;...
                eq.emp_share_6_10_by_skill_for;...
                eq.emp_share_11p_by_skill_for;...
                eq.emp_share_6_10_by_skill_inf;...
                eq.emp_share_11p_by_skill_inf];
            predictedNationalContributions = [eq.labor_share;inf_elast;...
                eq.share_formal_workers_500p/eq.share_formal_workers_100p];
        end        
        
        function numParametersSet = setModelParameters(...
                this,iRegion,x,paramNames)
            NP = length(paramNames);
            position = 0;
            for i = 1:NP
                position = position + 1;
                switch paramNames{i}
                    case {'alpha','gamma','eta',...
                            'sigmaFor','sigmaInf'}
                    %Parameters must lie inside (0,1)
                    v = 1/(1+exp(-x(position)));
                case {'xi_cons','xi_s','xi_rel_s','xi_rel_inf','A','coefRho1_sigma','coefRho2','B_z_ref','B_exp'}
                    %Parameters in (0,infty)
                    v = exp(x(position));
                case {'T';'rhoGridMax2Min'}
                    %Parameters in (1,infty)
                    v = exp(x(position)) + 1;
                case {'coefRho1_mu'}
                    %Unrestricted parameters
                    v = x(position);
                case {'ut_unemp'}
                    %2x1 vector of unrestricted parameters
                    v = x(position:(position+1));
                    position = position + 1;
                otherwise
                    error(['Unrecognized property/parameter/variable: ' ...
                        paramNames{i}]);
                end
                this.regionModels{iRegion}.(paramNames{i}) = v;
            end   
            numParametersSet = position;
        end
        
        function setAllParameters(this,x)
            for iRegion = 1:this.numRegions
                numParamsSet = this.setModelParameters(...
                    iRegion,x,this.regionSpecificParamNames);
                x(1:numParamsSet) = [];
            end
            for iRegion = 1:this.numRegions
                this.setModelParameters(...
                    iRegion,x,this.commonParamNames);
            end
        end
            
        
        function x = getModelTransformedParameters(this,iRegion,paramNames)
            NP = length(paramNames);
            x = [];
            for i = 1:NP
                v = this.regionModels{iRegion}.(paramNames{i});
                switch paramNames{i}
                    case {'alpha','gamma','eta',...
                            'sigmaFor','sigmaInf'}
                    %Parameters must lie inside (0,1)
                    x = [x; log(v/(1-v))];
                case {'A','xi_cons','xi_s','xi_rel_s','xi_rel_inf','coefRho1_sigma','coefRho2','B_z_ref','B_exp'}
                    %Parameters in (0,infty)
                    x = [x; log(v)];
                case {'T';'rhoGridMax2Min'}
                    %Parameters in (1,infty)
                    x = [x; log(v-1)];
                case {'coefRho1_mu','ut_unemp'}
                    %Unrestricted parameters
                    x = [x; v];
                otherwise
                    error(['Unrecognized property/parameter/variable: ' ...
                        paramNames{i}]);
                end
            end            
        end        
        
        function x = getAllParameters(this)
            x = [];
            for iRegion = 1:this.numRegions
                x = [x; this.getModelTransformedParameters(...
                    iRegion,this.regionSpecificParamNames)]; %#ok<AGROW>
            end
            x = [x; this.getModelTransformedParameters(...
                1,this.commonParamNames)];
        end
        
        function [regionalMoments, nationalContribs, eqs] = ...
                evaluateAllRegions(this,x,referenceEqs)
            if nargin < 3
                referenceEqs = repmat({[]},this.numRegions,1);
            end
            this.setAllParameters(x);
            regionalMoments = cell(this.numRegions,1);
            nationalContribs = cell(this.numRegions,1);
            eqs = cell(this.numRegions,1);
            for iRegion = 1:this.numRegions
                [regionalMoments{iRegion}, nationalContribs{iRegion}, eqs{iRegion}]...
                    = this.lossVector(iRegion,referenceEqs{iRegion});
            end
        end
        
        function [l,nationalMoments] = lossVal(this,regionalMoments,nationalContribs,verbose,x)
            if nargin < 5
                verbose = false;
            end
            allMoments = [];
            nationalMoments = zeros(3,1);
            for iRegion = 1:this.numRegions
                allMoments = [allMoments;regionalMoments{iRegion}];
                nationalMoments = nationalMoments + ...
                    this.regionSizes(iRegion) * nationalContribs{iRegion};
            end
            allMoments = [allMoments; nationalMoments];
            lossVector = allMoments - this.targetMoments;
            if this.weightOption == ModelEstimator.optimalWeights
                W = inv(this.momentsCovMatrix);
                W(end,end) = 1/this.momentsCovMatrix(end-1,end-1);
            elseif this.weightOption == ModelEstimator.sizeWeights
                w = [];
                nat_formality = 0;
                for iR = 1:this.numRegions
                    shareBySkill = [this.regionModels{iR}.eta;...
                        1-this.regionModels{iR}.eta];
                    skillWeights = repmat(shareBySkill,10,1);
                    share_1_5_pl_6_10 = (1-this.dataSummary{iR}.informal).*(...
                        1-this.dataSummary{iR}.fs_11p_for) + ...
                        this.dataSummary{iR}.informal.*(...
                        1-this.dataSummary{iR}.fs_11p_inf);
                    share_1_5_pl_11p = (1-this.dataSummary{iR}.informal).*(...
                        1-this.dataSummary{iR}.fs_6_10_for) + ...
                        this.dataSummary{iR}.informal.*(...
                        1-this.dataSummary{iR}.fs_6_10_inf);
                    momentWeights = [ones(8,1);...
                        share_1_5_pl_6_10;...
                        share_1_5_pl_11p;...
                        repmat(1-this.dataSummary{iR}.informal,2,1).*[this.dataSummary{iR}.fs_6_10_for;this.dataSummary{iR}.fs_11p_for];...
                        repmat(this.dataSummary{iR}.informal,2,1).*[this.dataSummary{iR}.fs_6_10_inf;this.dataSummary{iR}.fs_11p_inf]];
                    w = [w; this.regionSizes(iR) * ...
                        skillWeights.*momentWeights];
                    nat_formality = nat_formality + this.regionSizes(iR)*...
                        sum(shareBySkill.*(1-this.dataSummary{iR}.informal));
                end
                w = [w;1;1;nat_formality];
                W = diag(w);
                %Normalize all loss values, except the wage ones that
                % are already measured in relative terms
                mask = [repmat([true(4,1);false(8,1);true(8,1)],...
                    this.numRegions,1);true(3,1)];
                lossVector(mask) = lossVector(mask)./this.targetMoments(mask);
            else
                error('Invalid weightOption.');
            end
            l = sqrt(lossVector'*W*lossVector);
            %Display
            if ~verbose
                return;
            end
            message = ['x = [' num2str(x',16) ']\n' ...
                this.parameterString() 'RESIDUALS:\n'];
            regionalLosses = reshape(lossVector(1:(end-3)),20,this.numRegions)';
            nationalLosses = lossVector((end-2):end);
            for iRegion = 1:this.numRegions                    
                message = [message this.regionNames{iRegion} ...
                    ': ' num2str(regionalLosses(iRegion,:)) '\n'];
            end
            message = [message ...
                'Country target errors: ' num2str(nationalLosses') ...
                '\nLoss = ' ...
                num2str(l)];
            this.printMessage(message);
        end
        
        function gr = gradient(this,x,l,regionalMoments,nationalContribs,eqs)
            numCommonParams = length(this.getModelTransformedParameters(...
                1,this.commonParamNames));
            numRegionSpecificParams = length(this.getModelTransformedParameters(...
                1,this.regionSpecificParamNames));
            totalNumParams = numCommonParams + ...
                this.numRegions*numRegionSpecificParams;
            gr = nan(totalNumParams,1);
            col = 0;
            for iRegion = 1:this.numRegions
                for iParam = 1:numRegionSpecificParams %Only need to 
                    %re-calculate equilibrium for one region.
                    col = col + 1;
                    xp = x;
                    xp(col) = xp(col) + this.finiteDiffDelta;
                    this.setAllParameters(xp);      
                    [altRegionMom, altNatContrib] = this.lossVector(...
                        iRegion,eqs{iRegion});
                    altRegionalMoments = regionalMoments;
                    altRegionalMoments{iRegion} = altRegionMom;
                    altNationalContribs = nationalContribs;
                    altNationalContribs{iRegion} = altNatContrib;
                    lp = this.lossVal(altRegionalMoments,altNationalContribs);
                    gr(col) = (lp-l)/this.finiteDiffDelta;
                end
            end
            for i = 1:numCommonParams %Re-calculate everything
                col = col + 1;
                xp = x;
                xp(col) = xp(col) + this.finiteDiffDelta;
                [altRegionalMoments, altNationalContribs] = ...
                    this.evaluateAllRegions(xp,eqs);
                lp = this.lossVal(altRegionalMoments,altNationalContribs);
                gr(col) = (lp-l)/this.finiteDiffDelta;
            end
            this.printMessage(['Jacobian: ' num2str(gr') '\n']);
        end
        
        function [loss, gradient,nationalMoments,eqs] = estimationObjFun(this,x,calcGradient)
            failed = false;
            loss = [];
            gradient = [];
            try
                [regionalMoments, nationalContribs,eqs] = ...
                    this.evaluateAllRegions(x);
                [loss,nationalMoments] = this.lossVal(regionalMoments,nationalContribs,true,x);
                if isnan(loss) || isinf(loss) || (imag(loss)~=0)
                    failed = true;
                elseif calcGradient
                    gradient = this.gradient(x,loss,...
                        regionalMoments,nationalContribs,eqs);
                end
            catch excp
                failed = true;
                if ~strcmp(excp.identifier,'ModelEstimator:eqFailed')
                    this.printMessage(['Error calculating ' ...
                        'estimationObjFun.\n'...
                        getReport(excp)]);
                end
            end       
            allNums = [loss;gradient(:)];
            if any(isinf(allNums) | isnan(allNums))
                failed = true;
            end
            if failed
                loss = inf;
                numCommonParams = length(this.getModelTransformedParameters(...
                    1,this.commonParamNames));
                numRegionSpecificParams = length(this.getModelTransformedParameters(...
                    1,this.regionSpecificParamNames));
                totalNumParams = numCommonParams + ...
                    this.numRegions*numRegionSpecificParams;
                gradient = 10000*ones(totalNumParams,1);
                this.printMessage(['ModelEstimator:'...
                    'estimationObjFun failed with parameters:\n'...
                    this.parameterString()]);
            else
                this.printMessage(['ModelEstimator:'...
                    'estimationObjFun finished successfully.']);
            end
        end
        
        function printMessage(this,message,flourish)
            if nargin < 3
                flourish = false;
            end
            if isempty(this.messageFile)
                return;
            end
            if flourish
                starLine = repmat('x',1,50);
                starBlock = repmat([starLine '\n'],1,3);
                message = [starBlock message '\n' starBlock];
            end
            if strcmp(this.messageFile,'screen')
                fprintf([message '\n']);
                if strcmp(get(0,'Diary'),'on')
                    diary(get(0,'DiaryFile')); %Flush diary
                end
            else
                f=fopen(this.messageFile,'a');
                fprintf(f,[message '\n']);
                fclose(f);
            end
        end
        
        function initializeMessageFile(this)
            if isempty(this.messageFile) ...
                || strcmp(this.messageFile,'screen')
                return;
            end
            f=fopen(this.messageFile,'w');
            fprintf(f,['Log file: ' this.messageFile '\n\n']);
            fclose(f);
        end
        
        function str = parameterString(this)
            str = 'PARAMETERS:\n Common across regions: ';
            for i = 1:length(this.commonParamNames)
                param = this.commonParamNames{i};
                str = [str param '=' ...
                    num2str((this.regionModels{1}.(param))')];
                if i < length(this.commonParamNames)
                    str = [str ', '];
                end
            end
            for iRegion = 1:this.numRegions
                str = [str '\n ' this.regionNames{iRegion} ': '];
                for i = 1:length(this.regionSpecificParamNames)
                    param = this.regionSpecificParamNames{i};
                    model = this.regionModels{iRegion};
                    str = [str param '=' num2str((model.(param))')];
                    if i < length(this.regionSpecificParamNames)
                        str = [str ', '];
                    end
                end
            end
            str = [str '\n'];
        end
    end
end
