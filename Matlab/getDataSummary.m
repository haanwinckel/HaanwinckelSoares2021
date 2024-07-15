 function [summary,T] = getDataSummary(filename, etaByGroup)
    T = readtable(filename);
    
    if nargin < 2
        estimateEta = true;
        etaByGroup = [];
    else
        estimateEta = false;
    end
    
    regionCodes = unique(T.region);
 
    NR = length(regionCodes);
    
    ini_input = getStartingPoint(T, regionCodes, estimateEta);
    objFun = @(x) lossFunctionAllRegions(x,T,regionCodes, ...
        estimateEta,etaByGroup);
    opt = optimset('display','iter','maxiter',1000,'maxfunevals',10000);
    N_attempts = 5;
    minLoss = inf;
    for iAt = 1:N_attempts
        if iAt == 1
            ini = ini_input;
        else
            ini = ini_input + randn(size(ini_input));
        end
        [x,fv,fv2,flag,output] = lsqnonlin(objFun,ini,...
            [],[],opt);
        if fv < minLoss
            minLoss = fv;
            opt_x = x;
        end
    end
    
    [~, summary] = objFun(opt_x);
    totalSize = 0;
    for iR = 1:NR
        summary{iR}.eta = sum(T(T.region==regionCodes(iR),:).group_share...
            .*summary{iR}.eta_by_educ);
        summary{iR}.wh = exp(summary{iR}.lwh);

        sumSizeShares = summary{iR}.fs_1_5_for + summary{iR}.fs_6_10_for + summary{iR}.fs_11p_for;
        summary{iR}.fs_1_5_for = summary{iR}.fs_1_5_for./sumSizeShares;
        summary{iR}.fs_6_10_for = summary{iR}.fs_6_10_for./sumSizeShares;
        summary{iR}.fs_11p_for = summary{iR}.fs_11p_for./sumSizeShares;
        
        sumSizeShares = summary{iR}.fs_1_5_inf + summary{iR}.fs_6_10_inf + summary{iR}.fs_11p_inf;
        summary{iR}.fs_1_5_inf = summary{iR}.fs_1_5_inf./sumSizeShares;
        summary{iR}.fs_6_10_inf = summary{iR}.fs_6_10_inf./sumSizeShares;
        summary{iR}.fs_11p_inf = summary{iR}.fs_11p_inf./sumSizeShares;
      
        summary{iR}.regionSize = sum(...
            T(T.region==regionCodes(iR),:).group_share_all);     
        totalSize = totalSize + summary{iR}.regionSize;
    end
    for iR = 1:NR
        summary{iR}.regionSize = summary{iR}.regionSize/totalSize;
    end
end

function struct = vector2struct(input_region)    
    struct.eta_by_educ = [0;to01(input_region(1:10));1];
    rel_lambda_skill = [toPositive(input_region(11));1];
    rel_lambda_informal = toPositive(input_region(12));
    base_lambda_region = toPositive(input_region(13));
    struct.lambda_for = rel_lambda_skill*base_lambda_region;
    struct.lambda_inf = struct.lambda_for*rel_lambda_informal;
    struct.lwh = input_region(14:15);
    struct.wp_formal = input_region(16:17);
    struct.fs_1_5_for = to01(input_region(18:19));
    struct.fs_6_10_for = to01(input_region(20:21));
    struct.fs_11p_for = to01(input_region(22:23));
    struct.fs_1_5_inf = to01(input_region(24:25));
    struct.fs_6_10_inf = to01(input_region(26:27));
    struct.fs_11p_inf = to01(input_region(28:29));
    struct.wp_fs_6_10 = input_region(30:31);
    struct.wp_fs_11p = input_region(32:33);
    struct.phi = to01(input_region(34:35));
    struct.thetaq = toPositive(input_region(36:37));
end

function cellWithInputs = general2regionInputs(input,NR)
    numCommonInputs = 12;
    commonInputs = input(1:numCommonInputs);
    num_specific_inputs = (length(input)-numCommonInputs)/NR;
    cellWithInputs = cell(NR,1);
    for i = 1:NR
        firstRowSpecific = (i-1)*num_specific_inputs + numCommonInputs+1;
        lastRowSpecific = i*num_specific_inputs + numCommonInputs;
        cellWithInputs{i} = [commonInputs; input(...
            firstRowSpecific:lastRowSpecific,1)];
    end
end

function [residual, m] = lossFunctionSingleRegion(input_region,T_region)
    m=vector2struct(input_region);
    
    %From instantaneous hazards to monthly transitions
    for iSk = 1:2
        imp_formality = 1./...
            (1-(1-1/m.phi(iSk))*m.lambda_for(iSk)/m.lambda_inf(iSk));
        m.informal(iSk,1) = 1-imp_formality;

        dt = min(0.001,...
            0.9/max([m.lambda_for(iSk);m.lambda_inf(iSk);m.thetaq(iSk)]));
        trans_mat_inst =[...
            (1-m.thetaq(iSk)*dt)  m.thetaq(iSk)*dt*m.phi(iSk)  m.thetaq(iSk)*dt*(1-m.phi(iSk));...
            m.lambda_for(iSk)*dt 1-m.lambda_for(iSk)*dt       0          ;...
            m.lambda_inf(iSk)*dt        0         1-m.lambda_inf(iSk)*dt];
        trans_mat = trans_mat_inst^(1/dt);
        m.imp_F2U(iSk,1) = trans_mat(2,1);
        m.imp_I2U(iSk,1) = trans_mat(3,1);
        m.imp_U2F(iSk,1) = trans_mat(1,2);
        m.imp_U2I(iSk,1) = trans_mat(1,3);
        m.imp_U2E(iSk,1) = m.imp_U2F(iSk,1) + m.imp_U2I(iSk,1);
    end
    
    %Shares of workers in selected samples, by age-educ group
    %These will define mixing weights for skills in each age-educ group
    mean_lambda = (1-m.informal).*m.lambda_for+m.informal.*m.lambda_inf;
    m.unemp = mean_lambda./(mean_lambda+m.thetaq);
    emp = 1-m.unemp;
    shareSk_employed = m.eta_by_educ*emp(1)./...
        (m.eta_by_educ*emp(1)+(1-m.eta_by_educ)*emp(2));
    shareSk_unemployed = m.eta_by_educ*m.unemp(1)./...
        (m.eta_by_educ*m.unemp(1)+(1-m.eta_by_educ)*m.unemp(2));
    shareSk_formal = m.eta_by_educ*emp(1)*(1-m.informal(1))./...
        (m.eta_by_educ*emp(1)*(1-m.informal(1))...
        +(1-m.eta_by_educ)*emp(2)*(1-m.informal(2)));
    shareSk_informal = m.eta_by_educ*emp(1)*m.informal(1)./...
        (m.eta_by_educ*emp(1)*m.informal(1)...
        +(1-m.eta_by_educ)*emp(2)*m.informal(2));
    weights_employed = [shareSk_employed (1-shareSk_employed)];
    weights_unemployed = [shareSk_unemployed (1-shareSk_unemployed)];
    weights_formal = [shareSk_formal (1-shareSk_formal)];
    weights_informal = [shareSk_informal (1-shareSk_informal)];
    
    tab_predicted_moments = [...
        weights_employed * m.lwh ...
        weights_employed * m.wp_formal ...
        weights_formal * m.fs_1_5_for ...
        weights_formal * m.fs_6_10_for ...
        weights_formal * m.fs_11p_for ...
        weights_informal * m.fs_1_5_inf ...
        weights_informal * m.fs_6_10_inf ...
        weights_informal * m.fs_11p_inf ...
        weights_employed * m.wp_fs_6_10 ...
        weights_employed * m.wp_fs_11p ...
        weights_employed * m.informal ...
        weights_unemployed * m.imp_U2E ...
        weights_formal * m.imp_F2U ...
        weights_informal * m.imp_I2U ];
    
    fs_1_5_for = (1-T_region.fs_6_10_for-T_region.fs_11p_for);
    fs_1_5_inf = (1-T_region.fs_6_10_inf-T_region.fs_11p_inf);
    tab_data = [T_region.lwh T_region.wp_formal ...
        fs_1_5_for T_region.fs_6_10_for T_region.fs_11p_for ... 
        fs_1_5_inf T_region.fs_6_10_inf T_region.fs_11p_inf ... 
        T_region.wp_fs_6_10 T_region.wp_fs_11p ... 
        T_region.informal T_region.trans_U2E ...
        T_region.trans_F2U T_region.trans_I2U ];
    vert_weights = sqrt(T_region.group_share_all);
    hor_weights = [ones(1,2)*sum(T_region.group_share_all) ...
        sum(T_region.group_share_all.*(1-T_region.informal).*fs_1_5_for) ...
        sum(T_region.group_share_all.*(1-T_region.informal).*T_region.fs_6_10_for) ...
        sum(T_region.group_share_all.*(1-T_region.informal).*T_region.fs_11p_for) ...
        sum(T_region.group_share_all.*T_region.informal.*fs_1_5_inf) ...
        sum(T_region.group_share_all.*T_region.informal.*T_region.fs_6_10_inf) ...
        sum(T_region.group_share_all.*T_region.informal.*T_region.fs_11p_inf) ...
        sum(T_region.group_share_all.*(...
        (1-T_region.informal).*(fs_1_5_for+T_region.fs_6_10_for) + ...
        T_region.informal.*(fs_1_5_inf+T_region.fs_6_10_inf))) ...
        sum(T_region.group_share_all.*(...
        (1-T_region.informal).*(fs_1_5_for+T_region.fs_11p_for) + ...
        T_region.informal.*(fs_1_5_inf+T_region.fs_11p_inf))) ...
        sum(T_region.group_share_all) *ones(1,4)]/sum(T_region.group_share_all);
    hor_weights = 1;
    weights = vert_weights*hor_weights;
    
    residual = (tab_predicted_moments - tab_data).*sqrt(weights);
    residual = residual(:);
end

function [residual, structs] = lossFunctionAllRegions(input,T,regionCodes,...
    estimateEta, etaByGroup)
    if ~estimateEta
        input = [from01(etaByGroup(2:11));input];
    end
    NR = length(regionCodes);
    inputRegions = general2regionInputs(input,NR);
    residual = [];
    structs = cell(NR,1);
    for i = 1:NR
        T_region = T(T.region==regionCodes(i),:);
        [res, m] = lossFunctionSingleRegion(inputRegions{i},T_region);
        residual = [residual; res];
        structs{i} = m;
    end
end

function x = getRegionSpecificStartingPoint(T_region)
    
    x = [fromPositive(0.1);... Region-specific factor in E2U transitions
        T_region.lwh([12;1]);...
        T_region.wp_formal([12;1]);... 
        from01(1-T_region.fs_6_10_for([12;1])-T_region.fs_11p_for([12;1]))
        from01(T_region.fs_6_10_for([12;1]));... 
        from01(T_region.fs_11p_for([12;1]));... 
        from01(1-T_region.fs_6_10_inf([12;1])-T_region.fs_11p_inf([12;1]))
        from01(T_region.fs_6_10_inf([12;1]));... 
        from01(T_region.fs_11p_inf([12;1]));... 
        T_region.wp_fs_6_10([12;1]);... 
        T_region.wp_fs_11p([12;1]);... 
        from01(T_region.informal([12;1]));... 
        fromPositive(T_region.trans_U2E([12;1])+0.01)];
end

function x = getStartingPoint(T, regionCodes, estimateEta)
    if estimateEta
        ini_eta_2 = 1/2;
        x = from01(ini_eta_2)*ones(10,1);
    else
        x = [];
    end
    ini_rel_trans2U_high_skill = 0.5;
    ini_rel_trans2U_informal = 4;
    x = [x;...
        fromPositive([ini_rel_trans2U_high_skill;ini_rel_trans2U_informal])];
    NR = length(regionCodes);
    for i = 1:NR
        T_region = T(T.region==regionCodes(i),:);
        x = [x; ...
            getRegionSpecificStartingPoint(T_region)];
    end
end

function residual = transitionsLossFunction(...
    input, F2U, I2U, U2E, informality)

    thetaq = toPositive(input(1));
    phi = to01(input(2));
    lambda_for = toPositive(input(3));
    lambda_inf = toPositive(input(4));
    
    imp_formality = 1./...
        (1-(1-1/phi)*lambda_for/lambda_inf);
    imp_informality = 1-imp_formality;
    
    dt = min(0.001,...
        0.9/max([lambda_for;lambda_inf;thetaq]));
    trans_mat_inst =[...
        (1-thetaq*dt)  thetaq*dt*phi  thetaq*dt*(1-phi);...
        lambda_for*dt 1-lambda_for*dt       0          ;...
        lambda_inf*dt        0         1-lambda_inf*dt];
    trans_mat = trans_mat_inst^(1/dt);
    imp_F2U = trans_mat(2,1);
    imp_I2U = trans_mat(3,1);
    imp_U2F = trans_mat(1,2);
    imp_U2I = trans_mat(1,3);
    imp_U2E = imp_U2F + imp_U2I;
    
    residual = [imp_F2U;imp_I2U;imp_U2E;imp_informality] ...
        - [F2U;I2U;U2E;informality];     
end

function [thetaq, phi, lambda_for, lambda_inf] = ...
    monthly2instant(F2U, I2U, U2E, informality)
    objFun = @(x) transitionsLossFunction(x,...
        F2U,I2U,U2E,informality);
    opt = optimset('display','none','maxiter',1000,'maxfunevals',10000);
    [opt_x,fv,fv2,flag,output] = lsqnonlin(objFun,...
        zeros(4,1),[],[],opt);
    thetaq = toPositive(opt_x(1));
    phi = to01(opt_x(2));
    lambda_for = toPositive(opt_x(3));
    lambda_inf = toPositive(opt_x(4));
end

function y = fromPositive(x)
    y = log(x);
end

function x = toPositive(y)
    x = exp(y);
end

function y = from01(x)
    y = log(x./(1-x));
end

function x = to01(y)
    x = 1./(1+exp(-y));
end


