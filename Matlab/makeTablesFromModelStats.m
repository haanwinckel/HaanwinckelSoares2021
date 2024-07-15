function [ table_fit, table_cf, table_pol ] = ...
    makeTablesFromModelStats( resultsTable , makeFigures)

    if nargin < 2
        makeFigures = false;
    end

    moments = {'informality';'impliedUnemployment';'logWage';'wp_formal';...
        'wp_fs_6_10';'wp_fs_11p';...
        'fs_6_10_for';'fs_11p_for';'fs_6_10_inf';'fs_11p_inf'};
    titles = {'Imp. unemp.';'Informality';'log wage';'Formal wage premium';...
        'Size wage premium, 6-10';'Size wage premium, 11+';...
        'Sh. formal 6-10';'Sh. formal 11+';'Sh. informal 6-10';'Sh. informal 11+'};
    years = {'2003';'2012'};
    
    rowLabels = [];
    for iM = 1:length(moments)
        rowLabels = [rowLabels; moments(iM); {'    - Skilled'};...
            {'    - Unskilled'}];
    end
    
    %Table with fit    
    collate = @(a,b,c) num2cell(reshape([a; b; c],3*length(a),1));
    [ nat_all_data, nat_sk_data, nat_unsk_data, ...
        regions_all_data, regions_gap_data ] = ...
        getInfo( resultsTable, moments, 'data2003');
    [ nat_all_model, nat_sk_model, nat_unsk_model, ...
        regions_all_model, regions_gap_model, region_sizes ] = ...
        getInfo( resultsTable, moments, 'model2003');
    table_fit = [{'Moments' 'Data' 'Model'}; ...
        rowLabels collate(nat_all_data,nat_sk_data,nat_unsk_data) ...
        collate(nat_all_model,nat_sk_model,nat_unsk_model)];
%     
    %Make figures with fit
    if makeFigures
        f = makeFigWithFit(regions_all_data, regions_all_model,...
            regions_gap_data,regions_gap_model, nat_all_model, region_sizes, titles);
        saveFigure(f, '../Output/fig_fit.pdf', 12, 12);
    end
    
    %Table with changes and counterfactuals    
    [ nat_all_data12, nat_sk_data12, nat_unsk_data12, ...
        regions_all_data12, regions_gap_data12 ] = ...
        getInfo( resultsTable, moments, 'data2012');
    table_cf = [{'Moments' 'Data'}; ...
        rowLabels collate(nat_all_data12-nat_all_data,...
        nat_sk_data12-nat_sk_data,nat_unsk_data12-nat_unsk_data)];
    cols = {'model2012','cf_2012_ex_payroll_benefits','cf_2012_ex_rho','cf_2012_ex_mw','cf_2012_ex_educ','cf_2012_ex_A','cf_2012_ex_B1','cf_2003_plus_payroll_benefits','cf_2003_plus_rho','cf_2003_plus_mw','cf_2003_plus_educ','cf_2003_plus_A','cf_2003_plus_SBTC','model2012_robustness_T','model2012_robustness_coefRho1_mu'};
    for c = 1:length(cols)        
        [ nat_all, nat_sk, nat_unsk, regions_all, regions_gap ] = ...
            getInfo( resultsTable, moments, cols{c});
        table_cf = [table_cf [cols(c); ...
            collate(nat_all-nat_all_model,...
            nat_sk-nat_sk_model,nat_unsk-nat_unsk_model)]];
        if strcmp(cols{c},'model2012') && makeFigures
            f = makeFigWithFit(regions_all_data12-regions_all_data, ...
                regions_all-regions_all_model,...
                regions_gap_data12-regions_gap_data,...
                regions_gap-regions_gap_model, nat_all-nat_all_model, ...
                region_sizes, titles);
            saveFigure(f, ['../Output/fig_' cols{c} '.pdf'], 12, 12);
        end
    end
  
    %Table with policy simulations  
    moments = {'informality';'impliedUnemployment';'logWage';'laborShare';'govtSurplus';'netOutput'};
    
    rowLabels = [];
    for iM = 1:length(moments)
        if iM <= 3
            rowLabels = [rowLabels; moments(iM); {'    - Skilled'};...
                {'    - Unskilled'}];
        else
            rowLabels = [rowLabels; moments(iM)];
        end
    end  
    [ nat_all_model12, nat_sk_model12, nat_unsk_model12] = ...
        getInfo( resultsTable, moments, 'model2012' );
    table_pol = [{'Moments'}; rowLabels];
    cols = {'cf_policy_tau_minus_1pp','cf_policy_tau_u_minus_1pp','cf_policy_tau_u_minus_10pp','cf_policy_abono_duplo'};
    for c = 1:length(cols)        
        [ nat_all, nat_sk, nat_unsk] = ...
            getInfo( resultsTable, moments, cols{c});
        table_pol = [table_pol [cols(c); ...
            collate(nat_all(1:3)-nat_all_model12(1:3),...
            nat_sk(1:3)-nat_sk_model12(1:3),nat_unsk(1:3)-nat_unsk_model12(1:3));...
            nat_all(4)-nat_all_model12(4);...
            num2cell(nat_all(5:6)'./nat_all_model12(5:6)'-1)]];
    end
        
end

function [ nat_all, nat_sk, nat_unsk, regions_all, regions_gap, region_sizes ] = ...
    getInfo( resultsTable, moments, what)
    
    T = resultsTable(strcmp(resultsTable.what, what),:);
    
    %Get weights to aggregate measures at the region and then country
    %levels
    weights_sk = T(strcmp(T.skillLevel,'Skilled'),:).popShare;
    weights_unsk = T(strcmp(T.skillLevel,'Unskilled'),:).popShare;
    shSk_pop = [weights_sk weights_unsk]./(weights_sk+weights_unsk);
    unemp_skilled =  T(strcmp(T.skillLevel,'Skilled'),:).impliedUnemployment;
    unemp_unskilled =  T(strcmp(T.skillLevel,'Unskilled'),:).impliedUnemployment;
    shSk_emp = [weights_sk.*(1-unemp_skilled) weights_unsk.*(1-unemp_unskilled)]...
        ./(weights_sk.*(1-unemp_skilled)+weights_unsk.*(1-unemp_unskilled));
    shSk_unemp = [weights_sk.*unemp_skilled weights_unsk.*unemp_unskilled]...
        ./(weights_sk.*unemp_skilled+weights_unsk.*unemp_unskilled);
    inf_skilled =  T(strcmp(T.skillLevel,'Skilled'),:).informality;
    inf_unskilled =  T(strcmp(T.skillLevel,'Unskilled'),:).informality;
    shSk_formal = shSk_emp.*[(1-inf_skilled) (1-inf_unskilled)]...
        ./sum(shSk_emp.*[(1-inf_skilled) (1-inf_unskilled)],2);
    shSk_informal = shSk_emp.*[inf_skilled inf_unskilled]...
        ./sum(shSk_emp.*[inf_skilled inf_unskilled],2);
    for iM = 1:length(moments)
        regions_sk(:,iM) = T(strcmp(T.skillLevel,'Skilled'),:).(moments{iM});
        regions_unsk(:,iM) = T(strcmp(T.skillLevel,'Unskilled'),:).(moments{iM});
        switch moments{iM}
            case {'impliedUnemployment'}
                w = shSk_pop;
            case {'informality';'logWage';'wp_formal';'wp_fs_6_10';'wp_fs_11p'}
                w = shSk_emp;
            case {'fs_6_10_for';'fs_11p_for'}
                w = shSk_formal;
            case {'fs_6_10_inf';'fs_11p_inf'}
                w = shSk_informal;
        end
        regions_all(:,iM) = w(:,1).*regions_sk(:,iM)+w(:,2).*regions_unsk(:,iM);
    end
    regions_gap = regions_sk - regions_unsk;
    region_sizes = (weights_sk+weights_unsk)/sum(weights_sk+weights_unsk);
    nat_sk = sum(regions_sk.*region_sizes);
    nat_unsk = sum(regions_unsk.*region_sizes);
    nat_all = sum(regions_all.*region_sizes);
end

function f = makeFigWithFit(regions_all_data, regions_all_model,...
    regions_gap_data,regions_gap_model, nat_all_model, region_sizes, titles)
    %Make figures with fit
    f = figure;
    for iM = 1:10
        subplot(5,4,2*(iM-1)+1);
        vvv_vert = [min([regions_all_data(:,iM);regions_all_model(:,iM)]);
            max([regions_all_data(:,iM);regions_all_model(:,iM)])];
        l = diff(vvv_vert);
        vvv_vert = [vvv_vert(1)-0.1*l;vvv_vert(2)+0.1*l];
        vvv_hor = [min(regions_all_data(:,iM));
            max(regions_all_data(:,iM))];
        l = diff(vvv_hor);
        vvv_hor = [vvv_hor(1)-0.1*l;vvv_hor(2)+0.1*l];
        plot(vvv_vert,vvv_vert,'k--');
        hold on;
        plot(vvv_vert,[0;0],'Color',0.5*ones(3,1));
        %Trendline
        mdl = fitlm(regions_all_data(:,iM),regions_all_model(:,iM));
        plot(vvv_vert,mdl.predict(vvv_vert),'k:');
        scatter(regions_all_data(:,iM), regions_all_model(:,iM),...
            200*region_sizes,'black');
        xlabel('Data');
        ylabel('Model');
        axis([vvv_hor' vvv_vert']);
        title([titles{iM} ', All']);
        
        subplot(5,4,2*(iM-1)+2);
        vvv_vert = [min([regions_gap_data(:,iM);regions_gap_model(:,iM)]);
            max([regions_gap_data(:,iM);regions_gap_model(:,iM)])];
        l = diff(vvv_vert);
        vvv_vert = [vvv_vert(1)-0.1*l;vvv_vert(2)+0.1*l];
        vvv_hor = [min(regions_gap_data(:,iM));
            max(regions_gap_data(:,iM))];
        l = diff(vvv_hor);
        vvv_hor = [vvv_hor(1)-0.1*l;vvv_hor(2)+0.1*l];
        plot(vvv_vert,vvv_vert,'k--');        
        hold on;
        plot(vvv_vert,[0;0],...
            'Color',0.5*ones(3,1));
        %Trendline
        mdl = fitlm(regions_gap_data(:,iM),regions_gap_model(:,iM));
        plot(vvv_vert,mdl.predict(vvv_vert),'k:');
        scatter(regions_gap_data(:,iM), regions_gap_model(:,iM),...
            200*region_sizes,'black');
        xlabel('Data');
        ylabel('Model');
        axis([vvv_hor' vvv_vert']);
        title([titles{iM} ', sk. gap']);
    end
end
