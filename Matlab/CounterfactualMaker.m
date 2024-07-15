classdef CounterfactualMaker < handle
    
    properties (SetAccess = immutable, GetAccess = public)
        baselineModels;
        baselineEquilibria;        
        numRegions;
        regionNames;
        regionSizes;
        counterfactualList = containers.Map('KeyType','char','ValueType','any');
        y2012Models;
        y2012Equilibria;
        bootstrapId
    end
    properties (SetAccess = private, GetAccess = public)
        statTable = [];
        paramTable = [];
        y2012Calculated = false;
        currentModels;
        currentEquilibria;
        numericalOptions_fminunc = optimoptions('fsolve',...
            'display','iter',...
            'tolx',1e-3,'tolfun',1e-3,'maxiter',10000,'maxfunevals',10000);        
	numericalOptions_solveForChanges = optimoptions('fsolve',...
            'display','iter',...
            'tolx',1e-3,'tolfun',1e-3,'maxiter',20,'maxfunevals',10000);       
	maxAttempts_solveForChanges = 5;

    end
    properties (SetAccess = public, GetAccess = public)   
        messageFile = 'screen';
    end
    
    methods
        function this = CounterfactualMaker(estimator,estResults,...
                dataFile2012,mwChange,rhoChange,tau2012,bD2012,...
                bF2012,bV2012,bootstrapId, runRobustnessNoEduc)
            if nargin < 10
                this.bootstrapId = 0; %Point estimate
            else
                this.bootstrapId = bootstrapId;
            end
            if nargin < 11
                runRobustnessNoEduc = false;
            end
            this.baselineModels = estimator.regionModels;
            this.baselineEquilibria = estResults.eqs;
            this.numRegions = estimator.numRegions;
            this.regionNames = estimator.regionNames;
            for i = 1:this.numRegions
                this.baselineModels{i}.opt_equilibrium = ...
                    optimoptions('fsolve',...
                    'tolx',1e-8,'tolfun',1e-8,'display','none',...
                    'SpecifyObjectiveGradient',true,...
                    'MaxIter',100);                
                if this.bootstrapId > 0
                    this.baselineModels{i}.max_error_equilibrium = 1e-2;
                else
                    this.baselineModels{i}.max_error_equilibrium = 1e-6;
                end
                this.currentModels{i,1} = this.baselineModels{i}.clone();
                this.currentEquilibria{i,1} = this.baselineEquilibria{i};
                this.regionSizes(i,1) = estimator.dataSummary{i}.regionSize;
            end
            this.makeParameterTable([...
                estimator.commonParamNames;...
                estimator.regionSpecificParamNames;...
                {'lambda_for';'lambda_inf'}]);
            this.addCurrentModelResultsToTable('model2003');
            
            moments2012 = getDataSummary(dataFile2012,...
                estimator.dataSummary{1}.eta_by_educ);
            F=fieldnames(moments2012{1});
            for i=1:this.numRegions
                for m=1:length(F)
                    moments2012{i}.(F{m})=real(moments2012{i}.(F{m}));
                end
            end
            this.addDataToTable(estimator.dataSummary,'data2003');
            this.addDataToTable(moments2012,'data2012');
            
            %Create precursor counterfactual -- 2012 except productivity
            precursorName = 'cf_2012_ex_AB';
            this.addCounterfactualToList(precursorName, ...
                {'mw';'tau';'bD';'bF';'bV'}, ...
                {mwChange;tau2012;bD2012;bF2012;bV2012}, ...
                [1;0;0;0;0;0],[0;0;0;0;0;0],2003);
            this.setModelParameters(precursorName);
            %Add new share skilled
            for i=1:this.numRegions
                this.currentModels{i}.eta = moments2012{i}.eta;
                this.currentModels{i}.coefRho1_mu = ...
                    this.currentModels{i}.coefRho1_mu + rhoChange;
            end
            this.solveModels();
            this.addCurrentModelResultsToTable(precursorName);
            
            %Get change in A, B and store 2012 results
            for i = 1:this.numRegions
                meanLogWageData2003_bySkill = estimator.dataSummary{i}.lwh;
                meanLogWageData2012_bySkill = moments2012{i}.lwh;
                meanLogWageModel2003_bySkill = this.baselineEquilibria{i}.mean_lw;
                meanLogWageTargets{i,1} = meanLogWageModel2003_bySkill + ...
                    (meanLogWageData2012_bySkill-meanLogWageData2003_bySkill);
            end            
            this.solveForChanges(meanLogWageTargets);
            tableOnes = ones(height(this.paramTable)/this.numRegions,1);
            new_A = [];
            new_B = [];
            for i = 1:this.numRegions
                this.y2012Models{i,1} = this.currentModels{i}.clone();
                this.y2012Equilibria{i,1} = this.currentEquilibria{i,1};
                new_A = [new_A; this.currentModels{i}.A*tableOnes];
                new_B = [new_B; this.currentModels{i}.B_z_ref*tableOnes];
            end
            this.paramTable.new_A = new_A;
            this.paramTable.new_B_z_ref = new_B;
              
            if this.bootstrapId == 0    
                %Counterfactuals that are not bootstrapped:
                this.addCounterfactualToList('cf_2003_plus_educ', ...
                    {'eta'}, {1}, [0],[1],2003);
                this.addCounterfactualToList('cf_2003_plus_payroll_benefits', ...
                    {'tau';'bD';'bF';'bV'}, {1;1;1;1}, [0;0;0;0],[1;1;1;1],2003);
                this.addCounterfactualToList('cf_2003_plus_rho', ...
                    {'coefRho1_mu'}, {1}, [0],[1],2003);
                this.addCounterfactualToList('cf_2003_plus_A', ...
                    {'A'}, {1}, [0],[1],2003);
                this.addCounterfactualToList('cf_2003_plus_SBTC', ...
                    {'B_z_ref'}, {1}, [0],[1],2003);
                this.addCounterfactualToList('cf_2003_plus_mw', ...
                    {'mw'}, {1}, [0],[1],2003);
            end
            
            %Other CF            
            this.addCounterfactualToList('cf_2012_ex_payroll_benefits', ...
                {'tau';'bD';'bF';'bV'}, {1;1;1;1}, [1;1;1;1],[0;0;0;0],2012);
            this.addCounterfactualToList('cf_2012_ex_rho', ...
                {'coefRho1_mu'}, {1}, [1],[0],2012);
            this.addCounterfactualToList('cf_2012_ex_A', ...
                {'A'}, {1}, [1],[0],2012);
            this.addCounterfactualToList('cf_2012_ex_mw', ...
                {'mw'}, {1}, [1],[0],2012);
            this.addCounterfactualToList('cf_2012_ex_B1', ...
                {'B_z_ref'}, {1}, [1],[0],2012);
            %Policy
            tau = this.y2012Models{1}.tau;
            this.addCounterfactualToList('cf_policy_tau_minus_1pp', ...
                {'tau'}, {(tau-0.01)./tau}, [0],[1],2012);
            this.addCounterfactualToList('cf_policy_tau_u_minus_1pp', ...
                {'tau'}, {[1;(tau(2)-0.01)/tau(2)]}, [0],[1],2012);
            this.addCounterfactualToList('cf_policy_tau_u_minus_10pp', ...
                {'tau'}, {[1;(tau(2)-0.1)/tau(2)]}, [0],[1],2012);
            this.addCounterfactualToList('cf_policy_abono_duplo', ...
                {'bF'}, {[1;2]}, [0],[1],2012);
            
            this.runAllCounterfactuals();
            st = this.statTable;
            pt = this.paramTable;
            save(['tmp_cf_final_' num2str(this.bootstrapId) '.mat'],'st','pt')
            %Create 2012 except education
            precursorName = 'cf_2012_ex_educ';
            this.addCounterfactualToList(precursorName, ...
                {'eta'}, {1},1,0,2012);
            this.setModelParameters(precursorName);
            this.solveModels();
            this.addCurrentModelResultsToTable(precursorName);  
            if this.bootstrapId == 0 && runRobustnessNoEduc
                %Robustness with no educ: get targets
                for i = 1:this.numRegions                    
                    infTarget{i,1} = this.y2012Equilibria{i}.informality;
                    m2012noEduc{i} = this.currentModels{i}.clone();
                    eq2012noEduc{i} = this.currentEquilibria{i};
                end
                %Do robustness
                robParams = {'T';'coefRho1_mu'};
                for r = 1:2
                    for i = 1:this.numRegions
                        this.currentModels{i} = m2012noEduc{i}.clone();
                        this.currentEquilibria{i} = eq2012noEduc{i};
                    end
                    this.solveForChanges(meanLogWageTargets,infTarget,...
                        robParams{r});
                    new_A = [];
                    new_B = [];
                    new_other = [];
                    for i = 1:this.numRegions
                        new_A = [new_A; this.currentModels{i}.A*tableOnes];
                        new_B = [new_B; this.currentModels{i}.B_z_ref*tableOnes];
                        new_other = [new_other; ...
                            this.currentModels{i}.(robParams{r})*tableOnes];
                    end
                    this.paramTable.(['rob' num2str(r) '_A']) = new_A;
                    this.paramTable.(['rob' num2str(r) '_B']) = new_B;
                    this.paramTable.(['rob' num2str(r) '_' robParams{r}]) = new_other;
                end
            end
        end
        
        function statTable = runAllCounterfactuals(this)
            cfNames = keys(this.counterfactualList);
            for i = 1:length(cfNames)
                if ~strcmp(cfNames{i},'cf_2012_ex_AB') && ...
                        ~strcmp(cfNames{i},'cf_2012_ex_educ')
                    this.runSingleCounterfactual(cfNames{i});
                end
            end
            statTable = this.statTable;
        end
        
        function addCounterfactualToList(this,name,paramList,paramVals,...
                valRelativeTo2003,valRelativeTo2012,baseYear)
            %Unless otherwise specified, will use 2003 values
            cf.paramList = paramList;
            cf.paramVals = paramVals;
            cf.valRelativeTo2003 = valRelativeTo2003;            
            cf.valRelativeTo2012 = valRelativeTo2012;
            cf.baseYear = baseYear;
            this.counterfactualList(name) = cf;
        end        
        
        function setModelParameters(this, cfName)
            cf = this.counterfactualList(cfName);
            for r = 1:this.numRegions
                if cf.baseYear == 2003
                    this.currentModels{r} = this.baselineModels{r}.clone();
                    this.currentModels{r}.verbosity = 5;
                elseif cf.baseYear == 2012
                    this.currentModels{r} = this.y2012Models{r}.clone();
                else
                    error('Invalid base year.');
                end
                for i = 1:length(cf.paramList)
                    param = cf.paramList{i};
                    if cf.valRelativeTo2003(i)
                        ref = this.baselineModels{r}.(param);
                    elseif cf.valRelativeTo2012(i)
                        ref = this.y2012Models{r}.(param);
                    else
                        ref = 1;
                    end
                    this.currentModels{r}.(param) = ref .* ...
                        (cf.paramVals{i});
                end
            end                        
        end
        
        function solveModels(this)
            for r = 1:this.numRegions
                b_eq = this.baselineEquilibria{r};
                this.currentEquilibria{r} = ...
                    this.currentModels{r}.solveForEquilibrium(...
                    b_eq.rU,b_eq.q,b_eq);
                if isempty(this.currentEquilibria{r})
                    error(['Could not solve for equilibrium in ' ...
                        this.regionNames{r} '.']);
                end
            end
        end
        
        function addDataToTable(this,dataMoments,label)
            for i = 1:this.numRegions
                mom = dataMoments{i};
                bstrapId = repmat(this.bootstrapId,2,1);
                region = {this.regionNames{i};this.regionNames{i}};
                what = {label;label};
                skillLevel = {'Skilled';'Unskilled'};
                popShare = this.regionSizes(i) * ...
                    [mom.eta;1-mom.eta];
                logWage = mom.lwh;
                fs_6_10_for = mom.fs_6_10_for;
                fs_11p_for = mom.fs_11p_for;
                fs_6_10_inf = mom.fs_6_10_inf;
                fs_11p_inf = mom.fs_11p_inf;
                wp_fs_6_10 = mom.wp_fs_6_10;
                wp_fs_11p = mom.wp_fs_11p;
                wp_formal = mom.wp_formal;
                informality = mom.informal;
                instJobFindingRate = mom.thetaq;
                mean_lambda = (1-informality).*mom.lambda_for + ...
                    informality.*mom.lambda_inf;
                impliedUnemployment = mean_lambda ./ ...
                    (mean_lambda + instJobFindingRate);
                govtSurplus = nan(2,1);
                laborShare = nan(2,1);
                netOutput = nan(2,1);
                T = table(bstrapId,what,region,skillLevel,...
                    popShare,logWage,wp_formal,...
                    informality,instJobFindingRate,...                    
                    fs_6_10_for, fs_11p_for, fs_6_10_inf, fs_11p_inf, ...
                    wp_fs_6_10, wp_fs_11p, ...
                    impliedUnemployment,govtSurplus,laborShare,netOutput);
                this.statTable = [this.statTable; T];
            end
        end
        
        function addCurrentModelResultsToTable(this,label)
            for i = 1:this.numRegions
                eq = this.currentEquilibria{i};
                m = this.currentModels{i};
                bstrapId = repmat(this.bootstrapId,2,1);
                region = {this.regionNames{i};this.regionNames{i}};
                what = {label;label};
                skillLevel = {'Skilled';'Unskilled'};
                popShare = this.regionSizes(i) * ...
                    [m.eta;1-m.eta];
                logWage = eq.mean_lw;
                fs_6_10_for = eq.emp_share_6_10_by_skill_for;
                fs_11p_for = eq.emp_share_11p_by_skill_for;
                fs_6_10_inf = eq.emp_share_6_10_by_skill_inf;
                fs_11p_inf = eq.emp_share_11p_by_skill_inf;
                wp_fs_6_10 = eq.wp_fs_6_10;
                wp_fs_11p = eq.wp_fs_11p;
                wp_formal = eq.wp_formal;
                informality = eq.informality_by_skill;
                instJobFindingRate = eq.theta.*eq.q;
                impliedUnemployment = eq.unemp_by_skill;
                govtSurplus = repmat(eq.govt_surplus,2,1);
                laborShare = repmat(eq.labor_share,2,1);
                netOutput = repmat(eq.net_output,2,1);
                T = table(bstrapId,what,region,skillLevel,...
                    popShare,logWage,wp_formal,...
                    informality,instJobFindingRate,...                    
                    fs_6_10_for, fs_11p_for, fs_6_10_inf, fs_11p_inf, ...
                    wp_fs_6_10, wp_fs_11p, ...
                    impliedUnemployment,govtSurplus,laborShare,netOutput);
                this.statTable = [this.statTable; T];
            end
        end
        
         function makeParameterTable(this,paramNames)    
            for i = 1:this.numRegions
                m = this.currentModels{i};
                bstrapId = repmat(this.bootstrapId,2,1);
                region = {this.regionNames{i};this.regionNames{i}};
                skillLevel = {'Skilled';'Unskilled'};
                skillShare = [m.eta;1-m.eta];                
                T = table(bstrapId,region,skillLevel,skillShare);
                for p = 1:length(paramNames)
                    pn = paramNames{p};
                    val = m.(pn);
                    if length(val) == 1
                        val = [val;val];
                    end
                    if strcmp(pn,'coefRho1_mu')
                        T.coefRho1 = exp(val);
                    else
                        T.(pn) = val;
                    end
                end
                this.paramTable = [this.paramTable; T];
            end
        end
        
        function runSingleCounterfactual(this,cfName)
            disp(['Running counterfactual: ' cfName]);
            this.setModelParameters(cfName);
            this.solveModels();
            this.addCurrentModelResultsToTable(cfName);
        end
        
        function solveForChanges(this,meanLogWageTargets,informalityTarget, robustnessParam)
            %Assumes 2012 values have been added to this.currentModels
            if nargin < 3
                robustnessNoEducMode = false;
                robustnessParam = '';
            else
                robustnessNoEducMode = true;
                targets = [meanLogWageTargets; informalityTarget];
                ini_input = zeros(3,1);
            end
            
            for i = 1:this.numRegions
                baseParamVals{1} = this.baselineModels{i}.A;
                baseParamVals{2} = this.baselineModels{i}.B_z_ref;
                targets = meanLogWageTargets{i};
                ini_input = zeros(2,1);
                if robustnessNoEducMode
                    baseParamVals{3} = this.baselineModels{i}.(robustnessParam);
                    targets = [targets; informalityTarget{i}];
                    ini_input = zeros(3,1);
                end
                success = false;
                it = 0;
                best_fv = inf;
                while ~success && it < this.maxAttempts_solveForChanges
                    [xx,fv,flag,out] = fsolve(@(input) ...
                        this.solveForChangesObjFun(i,input,baseParamVals,targets,robustnessNoEducMode,robustnessParam),...
                        ini_input,this.numericalOptions_solveForChanges);
                    success = max(abs(fv)) < this.numericalOptions_solveForChanges.TolFun;
                    if success || max(abs(fv)) < best_fv
                        opt_inputs = xx;
                    end
                    it = it + 1;
                    ini_input = xx + 0.2*(rand(size(ini_input))-0.5);
                end
                this.solveForChangesObjFun(i,opt_inputs,baseParamVals,targets,robustnessNoEducMode,robustnessParam);
            end
            if robustnessNoEducMode
                this.addCurrentModelResultsToTable(['model2012_robustness_' robustnessParam]);
            else
                this.addCurrentModelResultsToTable('model2012');
            end                
        end
        
        function resid = solveForChangesObjFun(this,i,...
                inputChanges,baseParamVals,...
                targets, robustnessNoEducMode, robustnessParam)
            disp(['Entering solveForChangesObjFun, inputs: ' ...
                num2str(inputChanges')]);
            model = this.currentModels{i};
            v = model.verbosity;
            model.verbosity = 0;
            prevEq = this.currentEquilibria{i};
            model.A = baseParamVals{1} * exp(inputChanges(1));
            model.B_z_ref = baseParamVals{2} * exp(inputChanges(2)); 
            if robustnessNoEducMode
                model.(robustnessParam) = ...
                    baseParamVals{3} + inputChanges(3);
            end
            eq = model.solveForEquilibrium(prevEq.rU,...
                prevEq.q,prevEq);
            if isempty(eq)
                if robustnessNoEducMode
                    resid = inf(3,1);
                else
                    resid=inf(2,1);
                end
                return;
            end
            if robustnessNoEducMode
                vals=[eq.mean_lw;eq.informality];
            else
                vals = eq.mean_lw;
            end
            this.currentEquilibria{i} = eq;
            model.verbosity = v;
            
            resid = vals - targets;
            disp(['Exiting solveForChangesObjFun, resid: ' ...
                num2str(resid')]);
        end
    end
end
        
