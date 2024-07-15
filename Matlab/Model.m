classdef Model < handle
    
    properties
        %Note: all skill-specific parameters are 2x1 vectors
        % [skilled;unskilled]
        
        %Calibrated parameters
        eta = 0.5
        m = 0.0905
        tau = [0.7205;0.7205]
        bD = [7.48;4];        %Benefits: dismissal
        bF = [0.02;0.05];     %Benefits: fixed
        bV = [0.266;0.318];    %Benefits: variable
        lambda_for = [0.01;0.015];
        lambda_inf = [0.04;0.06];
        r = 0.008
        D = 0.30
        E = 0.50
        sigmaFor = 0.5;
        sigmaInf = 0.5;
        
        %Estimated parameters (the values below are defaults)
        A = 6
        B_z_ref = 0.25
        B_exp = 0.3
        alpha = 0.5
        gamma = 0.35
        coefRho1_mu = -1
        rhoGridMax2Min = 100000;
        coefRho2 = 1.75
        xi_cons = 4;
        xi_s = 0.1;
        xi_rel_s = 10;
        xi_rel_inf = 1;
        T = 2;
        ut_unemp = [-4;0];
        
        %Computed auxiliary parameters (for simplicity and efficiency)
        a
        b
        c
        tau_ratio
        sigma_tilde
        mw = 1
        
        %Computational properties
        zGrid;
        zGridSize = 10;
        max2minProd = 10000;
        rhoGridSize = 5;
        rhoGrid;
        rhoGridProbs;
        defaultIniLogNFor = @(z) repmat([-5+log(z)/2;0],1,3);
        defaultIniLogNInf = @(z) [-5+log(z)/2;0];
        useIniLogNFromPreviousIterInEq = true;
        previousIterEq = [];
        tolerance_integration = 1e-16
        opt_firm_problem = optimoptions('fsolve',...
            'tolx',1e-15,'tolfun',1e-15,'display','none');
        max_error_firm_problem = 1e-10;
        opt_equilibrium = optimoptions('fsolve',...
            'tolx',1e-10,'tolfun',1e-10,'display','none',...
            'SpecifyObjectiveGradient',true,...
            'MaxIter',100);
        opt_eq_base_alg = 'trust-region-dogleg';
        opt_eq_alt_alg = 'levenberg-marquardt';
        max_error_equilibrium = 1e-8;
        firmProblemMaxAttempts = 50;
        equilibriumMaxAttempts = 10;
        fpNumericalIntPoints = 1000;
        fpAux = [];
        
        verbosity = 0;
        outputLabel = ''; %This label is added to all output messages
        %Verbosity levels:
        % 0: No messages
        % 1: Errors only in solveForEquilibrium
        % 2: Details on solveForEquilibrium
        % 3: Errors in equilibriumConditions
        % 4: Details on equilibriumConditions
        % 5: Errors in firm problems
        % 6: Firm problem details
        verbosityOutput = 1; % A file handle, or 1 for output to screen        
    end
    
    methods
        
        function this = Model()
            this.calculate_auxiliary_params();
        end
        
        function calculate_auxiliary_params(this)
            this.a = 1+this.bV;
            this.b = (this.bF+this.bD.*this.lambda_for)*this.mw;
            this.c = (1-this.sigmaFor)*this.a+this.sigmaFor*(1+this.tau);
            this.tau_ratio = (1+this.tau)./this.a;
            this.sigma_tilde = (1+this.tau)*this.sigmaFor./this.c;
            this.zGrid = exp(linspace(0,log(this.max2minProd),this.zGridSize)');
            if this.rhoGridSize == 1
                this.rhoGrid = exp(this.coefRho1_mu);
                this.rhoGridProbs = 1;
            else
                logRhoVals = linspace(...
                    this.coefRho1_mu-log(this.rhoGridMax2Min)/2,...
                    this.coefRho1_mu+log(this.rhoGridMax2Min)/2,...
                    this.rhoGridSize)';
                this.rhoGrid = exp(logRhoVals);
                this.rhoGridProbs = ones(this.rhoGridSize,1)/this.rhoGridSize;
            end
            %fpAux is an auxiliary structure with pre-calculations
            % to speed up numerical integration in the problem of formal
            % firms.
            this.fpAux.eps = linspace(1e-14,1,this.fpNumericalIntPoints)';
            this.fpAux.eps_to_sigmaRatio_over_tr1 = ...
                this.fpAux.eps.^(((1-this.sigmaFor)/this.sigmaFor)...
                /this.tau_ratio(1));
            this.fpAux.eps_to_sigmaRatio_over_tr2 = ...
                this.fpAux.eps.^(((1-this.sigmaFor)/this.sigmaFor)...
                /this.tau_ratio(2));
            this.fpAux.eps_to_tr2_over_tr1 = ...
                this.fpAux.eps.^(this.tau_ratio(2)...
                /this.tau_ratio(1));
            this.fpAux.eps_to_tr1_over_tr2 = ...
                this.fpAux.eps.^(this.tau_ratio(1)...
                /this.tau_ratio(2));
            this.fpAux.spacing = diff(this.fpAux.eps(1:2));
        end
        
        function val = rho(this,ns,nu,rU,coefRho1)
            val = coefRho1*(ns+nu).^(this.coefRho2);
        end
        
        function val = rho_prime(this,ns,nu,rU,coefRho1)
            val = coefRho1*this.coefRho2*(ns+nu).^(this.coefRho2-1);
        end
        
        function val = B_fun(this,z)
            x = (z/this.B_z_ref).^this.B_exp;
            val = x./(1+x);
        end
        
        function val = z_cdf_unnormalized(this,z)
            val = 1 - (1./z).^this.T;
        end
        
        function val = z_density(this,z)
            norm = this.z_cdf_unnormalized(this.zGrid(end)) - ...
                this.z_cdf_unnormalized(this.zGrid(1));
            cons = this.T/norm;
            val = cons./z.^(this.T+1);
        end
        
        function val = xi_fun(this,v)
            val = this.xi_cons*(v + (exp(-this.xi_s*v)-1)/this.xi_s);
        end
        
        function val = xi_prime(this,v)
            val = this.xi_cons*(1-exp(-this.xi_s*v));
        end
        
        function val = upsilon_inf_Fi_rho(this,z,ns,nu,rU,coefRho1)
            ns_to_gamma = ns^this.gamma;
            nu_to_gamma = nu^this.gamma;
            B = this.B_fun(z);
            if this.gamma == 0
                F = this.A*z.*...
                    ns.^(this.alpha*B)...
                    .*nu.^(this.alpha*(1-B));
                Fs = this.alpha*B.*F./ns;
                Fu = this.alpha*(1-B).*F./nu;
            else
                bracket_term = (B.*ns_to_gamma...
                    +(1-B)*nu_to_gamma)...
                    ^(this.alpha/this.gamma-1);
                Fs = this.A*z.*B.*this.alpha...
                    *ns_to_gamma*bracket_term/ns;
                Fu = this.A*z.*(1-B)*this.alpha...
                    *nu_to_gamma*bracket_term/nu;
            end
            upsilon_F = 1/(1-this.sigmaInf*(1-this.alpha));
            upsilon_rho_rho_p = this.rho_prime(ns,nu,rU,coefRho1)...
                /(1+this.sigmaInf*(this.coefRho2-1));
            
            val = upsilon_F*[Fs;Fu]-upsilon_rho_rho_p;
        end
        
        function val = foc_inf(this,z,ns,nu,rU,q,coefRho1)
            rhs = rU + (this.r+this.lambda_inf) .* ...
                this.xi_prime(this.lambda_inf.*[ns;nu]./q)*this.xi_rel_inf.*[this.xi_rel_s;1] ...
                ./((1-this.sigmaInf)*q);
            val = this.upsilon_inf_Fi_rho(z,ns,nu,rU,coefRho1) - rhs;
        end
        
        function [n, w] = firm_problem_inf(this,grid_id,rU,q,ini_log_n)
            z = this.zGrid(grid_id);
            verbLevel_error = 5;
            verbLevel_detail = 6;
            tol = this.max_error_firm_problem;
            shock = 20;
            succeeded = true;
            for g = 1:this.rhoGridSize
                coefRho1 = this.rhoGrid(g);
                objFun = @(x) ...
                    this.foc_inf(z,exp(x(1)),exp(x(2)),rU,q,coefRho1);
                optim_ok = false;
                rep(g,1)=0;
                new_ini_log_n = ini_log_n{g}(z);
                while rep(g) < this.firmProblemMaxAttempts && ~optim_ok
                    ini_fv = objFun(new_ini_log_n);
                    ini_valid = ~any(isinf(ini_fv)|isnan(ini_fv));
                    if ini_valid
                        [log_n, fval, flag, output] = fsolve(objFun,...
                            new_ini_log_n,this.opt_firm_problem);
                        optim_ok = max(abs(fval)) < tol;                        
                        n{g,1} = exp(log_n);
                    end
                    rep(g,1) = rep(g,1)+1;
                    new_ini_log_n = ini_log_n{g}(z) + (rand(2,1)-0.5)*2*shock;
                end
                succeeded = succeeded && optim_ok;
            end
            
            if succeeded
                this.message(verbLevel_detail, ['firm_problem_inf, z = '...
                    num2str(z) ': solved (' ...
                    num2str(rep') ' attempts).']);
                for g = 1:this.rhoGridSize
                    coefRho1 = this.rhoGrid(g);
                    w{g,1} = this.bargained_wage_inf(z,n{g}(1),n{g}(2),rU,coefRho1);
                end
            else
                this.message(verbLevel_error, ['firm_problem_inf: ' ...
                    'could not solve for z=' num2str(z) ' (' ...
                    num2str(rep) ' attempts).']);
                for g = 1:this.rhoGridSize
                    n{g,1} = nan(2,1);
                    w{g,1} = nan(2,1);
                end
            end
        end
        
        function val = integrand_for_s_unc(this,z,ns,nu)
            B = this.B_fun(z);
            nsvec = this.fpAux.eps.*ns;
            ns_to_gamma = (nsvec).^(this.gamma);
            nu_to_gamma = (nu.*this.fpAux.eps_to_tr2_over_tr1).^(this.gamma);
            bracket_term = (B.*ns_to_gamma...
                +(1-B)*nu_to_gamma)...
                .^(this.alpha/this.gamma-1);
            Fs = this.A*z.*B.*this.alpha...
                *ns_to_gamma.*bracket_term./(nsvec);
            val = this.fpAux.eps_to_sigmaRatio_over_tr1.*Fs;
        end
        
        function val = integrand_for_u_unc(this,z,ns,nu)
            B = this.B_fun(z);
            nuvec = this.fpAux.eps.*nu;
            nu_to_gamma = (nuvec).^(this.gamma);
            ns_to_gamma = (ns.*this.fpAux.eps_to_tr1_over_tr2).^(this.gamma);
            bracket_term = (B.*ns_to_gamma...
                +(1-B)*nu_to_gamma)...
                .^(this.alpha/this.gamma-1);
            Fu = this.A*z.*(1-B)*this.alpha...
                *nu_to_gamma.*bracket_term./(nuvec);
            val = this.fpAux.eps_to_sigmaRatio_over_tr2.*Fu;
        end
        
        function val = upsilon_for_Fi_unc(this,z,ns,nu)
            val = [trapz(this.integrand_for_s_unc(...
                    z,ns,nu))*this.fpAux.spacing;
                trapz(this.integrand_for_u_unc(...
                    z,ns,nu))*this.fpAux.spacing;
                ] ./ this.sigma_tilde;
        end
        
        function val = foc_for_unc(this,z,ns,nu,rU,q)
            rhs = this.tau_ratio .* (rU-this.b) + ...
                (this.r+this.lambda_for).*...
                this.xi_prime(this.lambda_for.*[ns;nu]./q).*[this.xi_rel_s;1] ./ ...
                ((1-this.sigma_tilde).*q);
            val = this.upsilon_for_Fi_unc(z,ns,nu) - rhs;
        end
        
        function val = foc_for_strategic(this,z,ns,nu,rU,q)
            ups_Fi = this.upsilon_for_Fi_unc(z,ns,nu);
            rhs_s = this.tau_ratio(1) * (rU(1)-this.b(1)) + ...
                (this.r+this.lambda_for(1)).*this.xi_prime(...
                this.lambda_for(1).*ns./q(1))*this.xi_rel_s ./ ...
                ((1-this.sigma_tilde(1)).*q(1));
            foc_s = ups_Fi(1) - rhs_s;
            lhs_u = this.mw * this.c(2);
            rhs_u = (1-this.sigmaFor)*(rU(2)-this.b(2)) + ...
                this.sigmaFor * ups_Fi(2);
            foc_u = rhs_u - lhs_u;
            val = [foc_s;foc_u];
        end
        
        function val = integrand_for_strict(this,z,ns,nu)
            B = this.B_fun(z);
            nsvec = this.fpAux.eps.*ns;
            ns_to_gamma = (nsvec).^(this.gamma);
            nu_to_gamma = nu.^(this.gamma);
            Fs = this.A*z.*B.*this.alpha*ns_to_gamma.*...
                (B.*ns_to_gamma...
                +(1-B)*nu_to_gamma)...
                .^(this.alpha/this.gamma-1)...
                ./(nsvec);
            Fsu = Fs*(1-B)*(this.alpha-this.gamma)*...
                nu_to_gamma./(nu*...
                (B.*ns_to_gamma...
                +(1-B)*nu_to_gamma));   
            val = this.fpAux.eps_to_sigmaRatio_over_tr1.*[Fs Fsu];
        end
        
        function val = foc_for_strict(this,z,ns,nu,rU,q)
            B = this.B_fun(z);
            integral_terms = trapz(this.integrand_for_strict(...
                z,ns,nu))*this.fpAux.spacing;
            ups_Fs = integral_terms(1)/this.tau_ratio(1);
            rhs_s = this.sigma_tilde(1) * (rU(1)-this.b(1)) + ...
                this.sigmaFor*(this.r+this.lambda_for(1)).*...
                this.xi_prime(this.lambda_for(1).*ns./q(1))*this.xi_rel_s./ ...
                ((1-this.sigmaFor).*q(1));
            foc_s = ups_Fs - rhs_s;
            Fu = this.A*z.*(1-B)*this.alpha*...
                nu^(this.gamma-1) * ...
                (B.*ns^this.gamma...
                +(1-B)*nu^this.gamma)...
                .^(this.alpha/this.gamma-1);
            dws_dnu_plus_tax = integral_terms(2);
            rhs_u = this.mw*(1+this.tau(2)) + ...
                dws_dnu_plus_tax*ns + ...
                (this.r+this.lambda_for(2))*...
                this.xi_prime(this.lambda_for(2).*nu./q(2))./ ...
                q(2);
            foc_u = Fu - rhs_u;
            val = [foc_s;foc_u];
        end
        
        function val = bargained_wage_inf(this,z,ns,nu,rU,coefRho1)
            val = (1-this.sigmaInf)*rU + ...
                this.sigmaInf * this.upsilon_inf_Fi_rho(z,ns,nu,rU,coefRho1);
        end
        
        function val = bargained_wage_for_unc(this,z,ns,nu,rU)
            val = ((1-this.sigmaFor)*(rU-this.b) + ...
                this.sigmaFor * this.upsilon_for_Fi_unc(z,ns,nu))...
                ./ this.c;
        end
        
        function val = bargained_wage_s_for_strict(this,z,ns,nu,rU)
            B = this.B_fun(z);
            Fs = @(s,u) this.A*z.*B.*this.alpha*s.^(this.gamma-1).*...
                (B.*s.^this.gamma+(1-B)*u^this.gamma)...
                .^(this.alpha/this.gamma-1);
            integrand = @(eps) eps.^((1-this.sigmaFor)...
                /(this.sigmaFor*this.tau_ratio(1))).*Fs(eps*ns,nu);
            int_term = integral(integrand,1e-14,1,...
                'abstol',this.tolerance_integration,...
                'reltol',this.tolerance_integration);
            val = (1-this.sigmaFor)*(rU(1)-this.b(1))/this.c(1)...
                +int_term/(1+this.tau(1));
        end
        
        function [n, w] = ...
                firm_problem_for(this,grid_id,rU,q,ini_log_n_fun)
            z = this.zGrid(grid_id);
            ini_log_n = ini_log_n_fun(z);
            verbLevel_error = 5;
            verbLevel_detail = 6;
            tol = this.max_error_firm_problem;
            shock = 20;
            
            %Three possible kinds of solutions
            solutionTypeName = {'unconstrained';'strict';'strategic'};
            %With different objective functions....
            objFuns = {@(x) ...
                this.foc_for_unc(z,exp(x(1)),exp(x(2)),rU,q);...
                       @(x) ...
                this.foc_for_strict(z,exp(x(1)),exp(x(2)),rU,q);...
                       @(x) ...
                this.foc_for_strategic(z,exp(x(1)),exp(x(2)),rU,q)};
            %... and different wage functions
            wage_functions = {...
                @(n) this.bargained_wage_for_unc(z,n(1),n(2),rU);...
                @(n) [this.bargained_wage_s_for_strict(z,n(1),n(2),rU);...
                    this.mw];...
                @(n) this.bargained_wage_for_unc(z,n(1),n(2),rU)};
            %... and different starting points (columns in ini_log_n)
            
            %Solve
            rep = zeros(1,3);
            n = nan(2,3);
            w = nan(2,3);            
            for solutionType = 1:3
                lhs_strategic_u = this.mw * this.c(2);
                min_rhs_strategic_u = (1-this.sigmaFor)*(rU(2)-this.b(2));
                if (min_rhs_strategic_u >= lhs_strategic_u) && ...
                    ~strcmp('unconstrained',solutionTypeName{solutionType})
                    continue; %Outside option is so high that solution
                    % can never be strategic or binding
                end
                optim_ok = false;
                rep(solutionType)=0;
                new_ini_log_n = ini_log_n(:,solutionType);
                while rep(solutionType) < this.firmProblemMaxAttempts && ~optim_ok
                    ini_fv = objFuns{solutionType}(new_ini_log_n);
                    ini_valid = ~any(isinf(ini_fv)|isnan(ini_fv));
                    if ini_valid
                        [log_n, fval, flag, output] = ...
                            fsolve(objFuns{solutionType},...
                            new_ini_log_n,this.opt_firm_problem);
                        optim_ok = max(abs(fval)) < tol;
                    end
                    rep(solutionType) = rep(solutionType)+1;
                    new_ini_log_n = ini_log_n(:,solutionType)...
                        + (rand(2,1)-0.5)*2*shock;
                end
                if ~optim_ok
                    error(['firm_problem_for: ' ...
                        'could not solve for z=' num2str(z) ', '...
                        solutionTypeName{solutionType} '(' ...
                        num2str(rep(solutionType)) ' attempts).']);
                end
                n(:,solutionType) = exp(log_n);
                w(:,solutionType) = wage_functions{solutionType}(...
                    n(:,solutionType));
            end
            this.message(verbLevel_detail, ['firm_problem_for, z = '...
                num2str(z) ': solved (' ...
                num2str(rep) ' attempts).']);
        end
        
        function [w_for,w_inf] = wages(this,rU,q)
            w_for_unc = (1./this.a).*...
                (rU-this.b+(this.r+this.lambda_for)...
                .*this.sigmaFor.*this.xi./((1-this.sigmaFor)*q));
            w_for = [w_for_unc(1);max(w_for_unc(2),this.mw)];
            w_inf = rU + (this.r+this.lambda_inf)...
                .*this.sigmaInf.*this.xi./((1-this.sigmaInf)*q);
        end
        
        function val = F(this,z,ns,nu)
            B = this.B_fun(z);
            val = this.A*z.*(B.*ns.^(this.gamma)...
                +(1-B).*nu.^(this.gamma))...
                .^(this.alpha/this.gamma);
        end            
        
        function plotFormalProblem(this,z,rU,q,ini_ns,nu_min,nu_max)
            skillNames = {'Skilled';'Unskilled'};
            skillColors = {'k';'b'};
            focNames = {'Unconstrained','Constrained'};
            focStyles = {'-';'--'};
            legends = {};
            focFuns = {@(ns,nu) this.foc_for_unc(z,ns,nu,rU,q);...
                @(ns,nu) this.foc_for_strict(z,ns,nu,rU,q)};
            NP = 10;
            nu_vals = linspace(nu_min,nu_max,NP);
            figure; hold on;
            for iSkill = 1:2
                for iFocType = 1:2
                    fun = @(nu) fsolve(@(ns) sum(...
                        focFuns{iFocType}(ns,nu).*([1;2]==iSkill)),...
                        ini_ns,optimset('display','none'));
                    ns_vals = nan(NP,1);
                    for iP = 1:NP
                        ns_vals(iP) = fun(nu_vals(iP));
                    end
                    style = [skillColors{iSkill} focStyles{iFocType}];
                    plot(nu_vals,ns_vals,style);
                    legends = [legends; {[skillNames{iSkill} ' ' ...
                        focNames{iFocType}]}];
                end
            end
            %Min wage
            fun = @(nu) fsolve(@(ns) -this.mw + sum(...
                this.bargained_wage_for_unc(z,ns,nu,rU).*[0;1]),...
                ini_ns,optimset('display','none'));
            ns_vals = nan(NP,1);
            for iP = 1:NP
                ns_vals(iP) = fun(nu_vals(iP));
            end
            style = 'r';
            plot(nu_vals,ns_vals,style);
            legends = [legends; {'Minimum wage'}];
            legend(legends,'location''eastoutside');            
        end

        function fpsol = allFirmProblems(this,rU,q,...
                ini_log_n_for,ini_log_n_inf)
            verbLevel_error = 5;
            
            %Start by solving problem in the grid of z values
            I = length(this.zGrid);
            ngrid_for_s_all = nan(I,3);
            ngrid_for_u_all = nan(I,3);
            wgrid_for_s_all = nan(I,3);
            wgrid_for_u_all = nan(I,3);
            ngrid_inf_s = nan(I,this.rhoGridSize);
            ngrid_inf_u = nan(I,this.rhoGridSize);
            wgrid_inf_s = nan(I,this.rhoGridSize);
            wgrid_inf_u = nan(I,this.rhoGridSize);
            fpsol.formalSolutionType = nan(I,1);
            
            for i = 1:I
                [n, w] = this.firm_problem_for(i,rU,q,ini_log_n_for);
                ngrid_for_s_all(i,:) = n(1,:);
                ngrid_for_u_all(i,:) = n(2,:);
                wgrid_for_s_all(i,:) = w(1,:);
                wgrid_for_u_all(i,:) = w(2,:);
                [n, w] = this.firm_problem_inf(i,rU,q,ini_log_n_inf);
                for g = 1:this.rhoGridSize
                    ngrid_inf_s(i,g) = n{g}(1);
                    ngrid_inf_u(i,g) = n{g}(2);
                    wgrid_inf_s(i,g) = w{g}(1);
                    wgrid_inf_u(i,g) = w{g}(2);
                end
            end            
            
            %Create interpolants for employment choices
            method = 'pchip';
            %Informal is easy
            for g = 1:this.rhoGridSize
                pp_n_inf_s = interp1(log(this.zGrid),log(ngrid_inf_s(:,g)),method,'pp');
                pp_n_inf_u = interp1(log(this.zGrid),log(ngrid_inf_u(:,g)),method,'pp');
                pp_w_inf_s = interp1(log(this.zGrid),log(wgrid_inf_s(:,g)),method,'pp');
                pp_w_inf_u = interp1(log(this.zGrid),log(wgrid_inf_u(:,g)),method,'pp');
                fpsol.n_inf_s{g} = @(z) exp(ppval(pp_n_inf_s,log(z)));
                fpsol.n_inf_u{g} = @(z) exp(ppval(pp_n_inf_u,log(z)));
                fpsol.w_inf_s{g} = @(z) exp(ppval(pp_w_inf_s,log(z)));
                fpsol.w_inf_u{g} = @(z) exp(ppval(pp_w_inf_u,log(z)));
            end
            
            %For formal firms, you need to figure out which solution type
            %is the correct one.
            %First, an auxiliary function
            function v = getFormalVal(z,vals_interpolants)
                %Given a vector z of any size, will return a corresponding
                % vector of employment or wages from the valid regime.
                v = nan(size(z));
                strict_ids = z <= fpsol.formalHighestStrict;
                unc_ids = z >= fpsol.formalLowestUnc;
                strategic_ids = ~(strict_ids | unc_ids);
                v(unc_ids) = vals_interpolants{1}(z(unc_ids));
                v(strict_ids) = vals_interpolants{2}(z(strict_ids));
                v(strategic_ids) = vals_interpolants{3}(z(strategic_ids));
            end
            lhs_strategic_u = this.mw * this.c(2);
            min_rhs_strategic_u = (1-this.sigmaFor)*(rU(2)-this.b(2));
            if min_rhs_strategic_u >= lhs_strategic_u
                %rU_u is too high --- solution is surely unconstrained.
                fpsol.formalLowestUnc = [];
                fpsol.formalHighestStrict = [];
                pp_n_for_s = interp1(log(this.zGrid),log(ngrid_for_s_all(:,1)),method,'pp');
                pp_n_for_u = interp1(log(this.zGrid),log(ngrid_for_u_all(:,1)),method,'pp');
                pp_w_for_s = interp1(log(this.zGrid),log(wgrid_for_s_all(:,1)),method,'pp');
                pp_w_for_u = interp1(log(this.zGrid),log(wgrid_for_u_all(:,1)),method,'pp');

                fpsol.n_for_s = @(z) exp(ppval(pp_n_for_s,log(z)));
                fpsol.n_for_u = @(z) exp(ppval(pp_n_for_u,log(z)));
                fpsol.w_for_s = @(z) exp(ppval(pp_w_for_s,log(z)));
                fpsol.w_for_u = @(z) max(exp(ppval(pp_w_for_u,log(z))),this.mw);
                for tid = 1:3
                    fpsol.n_for_s_all{tid} = @(z) exp(ppval(pp_n_for_s,log(z)));
                    fpsol.n_for_u_all{tid} = @(z) exp(ppval(pp_n_for_u,log(z)));
                    fpsol.w_for_s_all{tid} = @(z) exp(ppval(pp_w_for_s,log(z)));
                    fpsol.w_for_u_all{tid} = @(z) max(exp(ppval(pp_w_for_u,log(z))),this.mw);
                end
            else
                %Interpolate all solution types
                for i = 1:3
                    pp_n_for_s = interp1(log(this.zGrid),log(ngrid_for_s_all(:,i)),method,'pp');
                    pp_n_for_u = interp1(log(this.zGrid),log(ngrid_for_u_all(:,i)),method,'pp');
                    pp_w_for_s = interp1(log(this.zGrid),log(wgrid_for_s_all(:,i)),method,'pp');
                    pp_w_for_u = interp1(log(this.zGrid),log(wgrid_for_u_all(:,i)),method,'pp');

                    fpsol.n_for_s_all{i} = @(z) exp(ppval(pp_n_for_s,log(z)));
                    fpsol.n_for_u_all{i} = @(z) exp(ppval(pp_n_for_u,log(z)));
                    fpsol.w_for_s_all{i} = @(z) exp(ppval(pp_w_for_s,log(z)));
                    fpsol.w_for_u_all{i} = @(z) max(exp(ppval(pp_w_for_u,log(z))),this.mw);
                end

                %Step 1: find lowest unconstrained firm
                leeway = @(z) [0 1]*this.bargained_wage_for_unc(z,...
                    fpsol.n_for_s_all{1}(real(z)),fpsol.n_for_u_all{1}(real(z)),rU)...
                    - this.mw;
                %Look for trouble:
                leewayGrid = arrayfun(leeway,this.zGrid);
                if any((leewayGrid(1:(end-1))>0) & (leewayGrid(2:end)<=0))
                    error(['allFirmProblems: disconnected segments in z ' ...
                        'space with unconstrained wages.']);
                end            
                lowIsUnc = leewayGrid(1) > 0;
                highIsUnc = leewayGrid(end) > 0;
                if lowIsUnc && highIsUnc %All unconstrained
                    fpsol.formalLowestUnc = this.zGrid(1);
                elseif highIsUnc %Some unconstrained
                    [fpsol.formalLowestUnc, fv] = binarySearch(...
                        leeway,this.zGrid(1),leewayGrid(1),...
                        this.zGrid(end),leewayGrid(end),...
                        this.opt_equilibrium.TolX,...
                        this.opt_equilibrium.TolX);
                    if abs(fv) > this.max_error_equilibrium
                        error(['allFirmProblems: Could not solve for '...
                            'lowest unconstrained formal firm. Residual: '...
                            num2str(fv)]);
                    end
                elseif lowIsUnc
                    error(['allFirmProblems: unexpected result, lowest ' ...
                        'is unconstrained but highest is not.']);
                else %None unconstrained
                    fpsol.formalLowestUnc = this.zGrid(end) + 1;
                end

                %Step 2: find highest strictly binding firm
                bindedness = @(z) this.mw - [0 1]*this.bargained_wage_for_unc(...
                    z,fpsol.n_for_s_all{2}(z),fpsol.n_for_u_all{2}(z),rU);
                %Look for trouble:
                bindednessGrid = arrayfun(bindedness,this.zGrid);
                crossings = bindednessGrid(1:(end-1)).*bindednessGrid(2:end)<0;
                if sum(crossings) > 1
                    error(['allFirmProblems: disconnected segments in z ' ...
                        'space with strictly binding wages.']);
                end            
                lowIsStrict = bindednessGrid(1) > 0;
                highIsStrict = bindednessGrid(end) > 0;
                if lowIsStrict && highIsStrict %All strictly binding
                    fpsol.formalHighestStrict = this.zGrid(end);
                elseif lowIsStrict %Some strictly binding
                    crossingLeftId = find(crossings);
                    [fpsol.formalHighestStrict, fv] = binarySearch(...
                        bindedness,this.zGrid(crossingLeftId),...
                        bindednessGrid(crossingLeftId),...
                        this.zGrid(crossingLeftId+1),...
                        bindednessGrid(crossingLeftId+1),...
                        this.opt_equilibrium.TolX,...
                        this.opt_equilibrium.TolX);
                    if abs(fv) > this.max_error_equilibrium
                        error(['allFirmProblems: Could not solve for '...
                            'highest strictly binding formal firm. ' ...
                            'Residual: ' num2str(fv)]);
                    end
                elseif highIsStrict
                    error(['allFirmProblems: unexpected result, highest is '...
                        'strictly binding but lowest is not.']);
                else %None strictly binding
                    fpsol.formalHighestStrict = this.zGrid(1) - 1;
                end

                if fpsol.formalHighestStrict > fpsol.formalLowestUnc
                    error(['allFirmProblems: unexpected result, highest '...
                        'strictly binding is above lowest unconstrained.']);
                end            
                fpsol.n_for_s = @(z) getFormalVal(z,fpsol.n_for_s_all);
                fpsol.n_for_u = @(z) getFormalVal(z,fpsol.n_for_u_all);
                fpsol.w_for_s = @(z) getFormalVal(z,fpsol.w_for_s_all);
                fpsol.w_for_u = @(z) getFormalVal(z,fpsol.w_for_u_all);
            end
            
            %Steady state instant profit values
            fpsol.pi_for = @ (z) ...
                this.F(z,fpsol.n_for_s(z),fpsol.n_for_u(z)) ...
                -sum([fpsol.n_for_s(z) fpsol.n_for_u(z)].*...
                (1+this.tau').*[fpsol.w_for_s(z) fpsol.w_for_u(z)]...
                +this.xi_fun([fpsol.n_for_s(z) fpsol.n_for_u(z)].*...
                this.lambda_for'./q').*[this.xi_rel_s 1],2);
            for i = 1:3                
                fpsol.pi_for_all{i,1} = @ (z) ...
                    this.F(z,fpsol.n_for_s_all{i}(z),fpsol.n_for_u_all{i}(z)) ...
                    -sum([fpsol.n_for_s_all{i}(z) fpsol.n_for_u_all{i}(z)].*...
                    (1+this.tau').*[fpsol.w_for_s_all{i}(z) fpsol.w_for_u_all{i}(z)]...
                    +this.xi_fun([fpsol.n_for_s_all{i}(z) fpsol.n_for_u_all{i}(z)].*...
                    this.lambda_for'./q').*[this.xi_rel_s 1],2);
            end
            for g = 1:this.rhoGridSize
                fpsol.pi_inf{g,1} = @ (z) ...
                    this.F(z,fpsol.n_inf_s{g}(z),fpsol.n_inf_u{g}(z)) ...
                    -this.rho(fpsol.n_inf_s{g}(z),fpsol.n_inf_u{g}(z),rU,this.rhoGrid(g))...
                    -sum([fpsol.n_inf_s{g}(z) fpsol.n_inf_u{g}(z)].*...
                    [fpsol.w_inf_s{g}(z) fpsol.w_inf_u{g}(z)] ...
                    +this.xi_fun([fpsol.n_inf_s{g}(z) fpsol.n_inf_u{g}(z)].*...
                    this.lambda_inf'./q').*[this.xi_rel_s 1]*this.xi_rel_inf,2);
            end            
            fpsol.informal_type_g = @(z,g) fpsol.pi_for(z)<fpsol.pi_inf{g}(z);
        end       
        
        function [rU_error,theta_error,eq] = equilibriumConditions(...
                this,rU,q,startingPointEq)            
            verbLevel_error = 3;
            verbLevel_details = 4;
            
            this.calculate_auxiliary_params();
            
            if nargin < 4
                startingPointEq = [];
            end
            [ini_log_n_for,ini_log_n_inf] = ...
                this.getInitialPoints(startingPointEq);
            
            this.message(verbLevel_details, ['eqConditions: ' ...
                'starting with rU=[' num2str(rU',16) '], q=[' ...
                num2str(q',16) '].']);
            warning('off','MATLAB:integral:MinStepSize');
            warning('off','MATLAB:integral:MaxIntervalCountReached');
            try
                fp = this.allFirmProblems(rU,q,...
                    ini_log_n_for,ini_log_n_inf);
            catch excp
                this.message(verbLevel_error, ['eqConditions: ' ...
                    'error when solving firm problems:\n' ...
                    getReport(excp)]);
                fp = [];
            end
            
            if isempty(fp)
                [rU_error,theta_error] = deal(inf(2,1));
                eq = [];
                this.message(verbLevel_error, ['eqConditions: ' ...
                    'Firm problem was not solved correctly.']);
                return;
            end
            
            tol_int = this.opt_firm_problem.TolX;
            
            N = zeros(4,1);%cols:for_s,for_u,inf_s,inf_u
            W = zeros(4,1);
            
            %Find firms that are indifferent between sectors
            lzv = linspace(log(this.zGrid(1)),log(this.zGrid(end)),100)';
            for g = 1:this.rhoGridSize
                profGapFun=@(lz) fp.pi_for(exp(lz))-fp.pi_inf{g}(exp(lz));
                profGap = profGapFun(lzv);
                id_indiff = find(profGap(1:(end-1)).*profGap(2:end)<0);
                fp.z_indiff{g} = zeros(length(id_indiff),1);
                for i = 1:length(id_indiff)
                    [sol,fv] = binarySearch(profGapFun,...
                        lzv(id_indiff(i)),profGap(id_indiff(i)),...
                        lzv(id_indiff(i)+1),profGap(id_indiff(i)+1),...
                        this.opt_equilibrium.TolX/1000,...
                        this.opt_equilibrium.TolX);
                    fp.z_indiff{g}(i,1) = exp(sol);
                end
                waypoints = sort([this.zGrid(1);this.zGrid(end);fp.z_indiff{g};...
                    fp.formalHighestStrict;fp.formalLowestUnc]);
                waypoints(waypoints<this.zGrid(1))=[];
                waypoints(waypoints>this.zGrid(end))=[];

                num_intervals = (length(waypoints)-1);
                for intId = 1:num_intervals
                    low_z = waypoints(intId);
                    high_z = waypoints(intId+1);
                    mid_z = (low_z+high_z)/2;
                    intOverFirms = @(fun) this.rhoGridProbs(g)*this.m*...
                        integral(fun,low_z,high_z,...
                        'abstol',tol_int,'reltol',tol_int);
                    if fp.informal_type_g(mid_z,g) %Informal firms in interval
                        N(3) = N(3) + intOverFirms(@(z) ( ...
                            fp.n_inf_s{g}(z').*this.z_density(z'))');
                        N(4) = N(4) + intOverFirms(@(z) ( ...
                            fp.n_inf_u{g}(z').*this.z_density(z'))');
                        W(3) = W(3) + intOverFirms(@(z) ( ...
                            fp.n_inf_s{g}(z').*fp.w_inf_s{g}(z').*this.z_density(z'))');
                        W(4) = W(4) + intOverFirms(@(z) (...
                            fp.n_inf_u{g}(z').*fp.w_inf_u{g}(z').*this.z_density(z'))');
                    else
                        if mid_z > fp.formalLowestUnc %Unconstrained
                            tid = 1;
                        elseif mid_z < fp.formalHighestStrict
                            tid = 2;
                        else
                            tid = 3;
                        end
                        N(1) = N(1) + intOverFirms(@(z) ( ...
                            fp.n_for_s_all{tid}(z').*this.z_density(z'))');
                        N(2) = N(2) + intOverFirms(@(z) ( ...
                            fp.n_for_u_all{tid}(z').*this.z_density(z'))');
                        W(1) = W(1) + intOverFirms(@(z) ( ...
                            fp.n_for_s_all{tid}(z').*fp.w_for_s_all{tid}(z').*this.z_density(z'))');
                        W(2) = W(2) + intOverFirms(@(z) ( ...
                            fp.n_for_u_all{tid}(z').*fp.w_for_u_all{tid}(z').*this.z_density(z'))');
                    end                    
                end  
            end

            eq.N_for_s = N(1);
            eq.N_for_u = N(2);
            eq.N_inf_s = N(3);
            eq.N_inf_u = N(4);

            eq.theta = (q/this.D).^(-1/this.E);
            theta_error = ...
                ((this.lambda_for+eq.theta.*q).*[eq.N_for_s;eq.N_for_u]...
                +(this.lambda_inf+eq.theta.*q).*[eq.N_inf_s;eq.N_inf_u])...
                ./(eq.theta.*q.*[this.eta;1-this.eta]) -1;
            
            eq.mean_w_for = W(1:2)./N(1:2);
            eq.mean_w_inf = W(3:4)./N(3:4);
            eq.mean_w_for(isnan(eq.mean_w_for))=1; %when phi == 0
            eq.mean_w_inf(isnan(eq.mean_w_inf))=1; %when phi == 1
            
            eq.E_for = (this.a.*eq.mean_w_for+this.b...
                +this.lambda_for.*rU/this.r)./(this.r+this.lambda_for);
            eq.E_inf = (eq.mean_w_inf...
                +this.lambda_inf.*rU/this.r)./(this.r+this.lambda_inf);
            
            eq.phi = [this.lambda_for(1)*eq.N_for_s...
                /(this.lambda_for(1)*eq.N_for_s+this.lambda_inf(1)*eq.N_inf_s);...
                this.lambda_for(2)*eq.N_for_u...
                /(this.lambda_for(2)*eq.N_for_u+this.lambda_inf(2)*eq.N_inf_u)];
            
            vals = [this.ut_unemp ...
                (this.a.*eq.mean_w_for+this.b) ...
                eq.mean_w_inf];
            weights = [1./(eq.theta.*q) ...
                eq.phi./(this.r+this.lambda_for) ...
                (1-eq.phi)./(this.r+this.lambda_inf)];
            imp_rU = sum(vals.*weights,2)./sum(weights,2);
            rU_error = imp_rU - rU;
            this.message(verbLevel_details, ['eqConditions: ' ...
                'finished with rU_error=[' num2str(rU_error',16) '], ' ...
                'theta_error=[' num2str(theta_error',16) '].']);
            
            eq.rU = rU;
            eq.q = q;
            eq.fp = fp;
            
            warning('on','MATLAB:integral:MinStepSize');
            warning('on','MATLAB:integral:MaxIntervalCountReached');
        end
        
        function eq = addEquilibriumStatistics(this,eq)
            warning('off','MATLAB:integral:MinStepSize');
            warning('off','MATLAB:integral:MaxIntervalCountReached');
            
            tol_int = this.opt_firm_problem.TolX;
            %The following variables are indexed as following:
            % i: skill
            % j: sector
            % k: firm size categories
            fsThresholds = [0;5.5;10.5;100.5;500.5;inf];
            nFS = length(fsThresholds)-1;
            [N, LW] = deal(zeros(2,2,nFS));
            PI = zeros(1,2); %At the sector level only
            fp=eq.fp;
            
            for g = 1:this.rhoGridSize
                waypoints = sort([this.zGrid(1);this.zGrid(end);fp.z_indiff{g};...
                    fp.formalHighestStrict;fp.formalLowestUnc]);
                waypoints(waypoints<this.zGrid(1))=[];
                waypoints(waypoints>this.zGrid(end))=[];

                num_intervals = (length(waypoints)-1);
                for intId = 1:num_intervals
                    low_z = waypoints(intId);
                    high_z = waypoints(intId+1);
                    mid_z = (low_z+high_z)/2;
                    if fp.informal_type_g(mid_z,g) %Informal firms in interval
                        j = 2;
                        n_s = fp.n_inf_s{g};
                        n_u = fp.n_inf_u{g};
                        lw_s = @(z) log(fp.w_inf_s{g}(z)).*fp.n_inf_s{g}(z);
                        lw_u = @(z) log(fp.w_inf_u{g}(z)).*fp.n_inf_u{g}(z);
                        pi = fp.pi_inf{g};
                    else
                        if mid_z > fp.formalLowestUnc %Unconstrained
                            tid = 1;
                        elseif mid_z < fp.formalHighestStrict
                            tid = 2;
                        else
                            tid = 3;
                        end
                        j = 1;
                        n_s = fp.n_for_s_all{tid};
                        n_u = fp.n_for_u_all{tid};
                        lw_s = @(z) log(fp.w_for_s_all{tid}(z)).*fp.n_for_s_all{tid}(z);
                        lw_u = @(z) log(fp.w_for_u_all{tid}(z)).*fp.n_for_u_all{tid}(z);
                        pi = fp.pi_for_all{tid};
                    end
                    %Check firm size categories
                    cat_highest_size = n_s(high_z)+n_u(high_z);
                    cat_lowest_size = n_s(low_z)+n_u(low_z);
                    for k = 1:nFS
                        if cat_lowest_size >= fsThresholds(k+1)
                            continue;
                        elseif cat_highest_size <= fsThresholds(k)
                            continue;
                        end
                        if cat_lowest_size < fsThresholds(k)
                            objFun = @(z) n_s(exp(z))+n_u(exp(z))-fsThresholds(k);
                            [lint_low_z, fv] = binarySearch(objFun,...
                                log(low_z),objFun(log(low_z)),...
                                log(high_z),objFun(log(high_z)),...
                                this.opt_equilibrium.TolX/10000,...
                                this.opt_equilibrium.TolX);
                            int_low_z = exp(lint_low_z);
                            if abs(fv)>this.max_error_equilibrium
                                error('Could not find integration limits.');
                            end
                        else
                            int_low_z = low_z;
                        end
                        if cat_highest_size > fsThresholds(k+1)
                            objFun = @(z) n_s(exp(z))+n_u(exp(z))-fsThresholds(k+1);
                            [lint_high_z, fv] = binarySearch(objFun,...
                                log(low_z),objFun(log(low_z)),...
                                log(high_z),objFun(log(high_z)),...
                                this.opt_equilibrium.TolX/10000,...
                                this.opt_equilibrium.TolX);
                            int_high_z = exp(lint_high_z);
                            if abs(fv)>this.max_error_equilibrium
                                error('Could not find integration limits.');
                            end
                        else
                            int_high_z = high_z;
                        end                            
                        intOverFirms = @(fun) this.rhoGridProbs(g)*this.m*...
                            integral(@(z) (fun(z').*this.z_density(z'))',...
                            int_low_z,int_high_z,...
                            'abstol',tol_int,'reltol',tol_int);
                        N(1,j,k) = N(1,j,k) + intOverFirms(n_s);
                        N(2,j,k) = N(2,j,k) + intOverFirms(n_u);
                        LW(1,j,k) = LW(1,j,k) + intOverFirms(lw_s);
                        LW(2,j,k) = LW(2,j,k) + intOverFirms(lw_u);
                        PI(j) = PI(j) + intOverFirms(pi);
                    end
                end  
            end            
            
            eq.unemp = 1-eq.N_for_s-eq.N_for_u-eq.N_inf_s-eq.N_inf_u;
            eq.unemp_by_skill = [this.eta-eq.N_for_s-eq.N_inf_s;...
                1-this.eta-eq.N_for_u-eq.N_inf_u]...
                ./[this.eta;1-this.eta];
            eq.informality = (eq.N_inf_s+eq.N_inf_u)/(1-eq.unemp);
            eq.informality_by_skill = [eq.N_inf_s;eq.N_inf_u]...
                ./([this.eta;1-this.eta].*(1-eq.unemp_by_skill));
            eq.total_profit_for = PI(1);
            eq.total_profit_inf = PI(2);
            eq.total_benefits = sum([eq.N_for_s;eq.N_for_u].*...
                ((this.a-1).*eq.mean_w_for+this.b));
            eq.total_taxes = sum([eq.N_for_s;eq.N_for_u].*this.tau...
                .*eq.mean_w_for);
            eq.govt_surplus = eq.total_taxes - eq.total_benefits;
            eq.net_output = eq.govt_surplus + ...
                eq.total_profit_for + eq.total_profit_inf + ...
                sum([eq.N_for_s;eq.N_for_u;eq.N_inf_s;eq.N_inf_u]...
                .*[eq.mean_w_for;eq.mean_w_inf]);
            eq.labor_share = ...
                sum([eq.N_for_s;eq.N_for_u;eq.N_inf_s;eq.N_inf_u]...
                .*[eq.mean_w_for;eq.mean_w_inf])...
                /(eq.net_output - eq.govt_surplus);
            %Firm size variables
            eq.emp_share_6_10_by_skill_for = [...
                N(1,1,2)/eq.N_for_s;...
                N(2,1,2)/eq.N_for_u];
            eq.emp_share_11p_by_skill_for = [...
                sum(N(1,1,3:end),3)/eq.N_for_s;...
                sum(N(2,1,3:end),3)/eq.N_for_u];
            eq.emp_share_6_10_by_skill_inf = [...
                N(1,2,2)/eq.N_inf_s;...
                N(2,2,2)/eq.N_inf_u];
            eq.emp_share_11p_by_skill_inf = [...
                sum(N(1,2,3:end),3)/eq.N_inf_s;...
                sum(N(2,2,3:end),3)/eq.N_inf_u];
            eq.emp_share_6_10_by_skill = [...
                (N(1,1,2)+N(1,2,2))/(eq.N_for_s+eq.N_inf_s);...
                (N(2,1,2)+N(2,2,2))/(eq.N_for_u+eq.N_inf_u)];
            eq.emp_share_11p_by_skill = [...
                sum(N(1,1,3:end)+N(1,2,3:end),3)/(eq.N_for_s+eq.N_inf_s);...
                sum(N(2,1,3:end)+N(2,2,3:end),3)/(eq.N_for_u+eq.N_inf_u)];
            eq.share_formal_workers_100p = ...
                sum(sum(N(:,1,4:end),1),3)/(eq.N_for_s+eq.N_for_u);
            eq.share_formal_workers_500p = ...
                sum(N(:,1,5),1)/(eq.N_for_s+eq.N_for_u);
            eq.mean_lw_for = sum(LW(:,1,:),3)./sum(N(:,1,:),3);
            eq.mean_lw_inf = sum(LW(:,2,:),3)./sum(N(:,2,:),3);
            eq.mean_lw = sum(sum(LW,3),2)./sum(sum(N,3),2);
            
            eq.wp_formal = exp(eq.mean_lw_for-eq.mean_lw_inf)-1;
            eq.wp_fs_6_10 = exp(...
                sum(LW(:,:,2),2)./sum(N(:,:,2),2)...
                -sum(LW(:,:,1),2)./sum(N(:,:,1),2))-1;
            eq.wp_fs_11p = exp(...
                sum(sum(LW(:,:,3:end),3),2)...
                ./sum(sum(N(:,:,3:end),3),2)...
                -sum(LW(:,:,1),2)./sum(N(:,:,1),2))-1;
            
            warning('on','MATLAB:integral:MinStepSize');
            warning('on','MATLAB:integral:MaxIntervalCountReached');
        end
        
        function [ini_log_n_for,ini_log_n_inf] = ...
                getInitialPoints(this,referenceEq)
            if isempty(referenceEq)
                ini_log_n_for = this.defaultIniLogNFor;
                for g = 1:this.rhoGridSize
                    ini_log_n_inf{g,1} = this.defaultIniLogNInf;
                end
            else
                if any(strcmp(fieldnames(referenceEq.fp),'n_for_s_all'))
                    ini_log_n_for = @(z) log([...
                        referenceEq.fp.n_for_s_all{1}(z) ...
                        referenceEq.fp.n_for_s_all{2}(z) ...
                        referenceEq.fp.n_for_s_all{3}(z);...
                        referenceEq.fp.n_for_u_all{1}(z) ...
                        referenceEq.fp.n_for_u_all{2}(z) ...
                        referenceEq.fp.n_for_u_all{3}(z)]);
                else
                    ini_log_n_for = @(z) log([...
                        referenceEq.fp.n_for_s(z) ...
                        referenceEq.fp.n_for_s(z) ...
                        referenceEq.fp.n_for_s(z);...
                        referenceEq.fp.n_for_u(z) ...
                        referenceEq.fp.n_for_u(z) ...
                        referenceEq.fp.n_for_u(z)]);
                end
                for g = 1:this.rhoGridSize
                    ini_log_n_inf{g,1} = @(z) log([referenceEq.fp.n_inf_s{g}(z);...
                        referenceEq.fp.n_inf_u{g}(z)]);
                end
            end
        end
        
        function [val,jac,rU,q] = equilibriumObjFun(this,x,...
                referenceEq)
            if nargin < 3 || isempty(referenceEq)
                if this.useIniLogNFromPreviousIterInEq ...
                        && ~isempty(this.previousIterEq)
                    referenceEq = this.previousIterEq;
                else
                    referenceEq = [];
                end
            end
            rU = exp(x(1:2));
            q = exp(x(3:4));
            [rU_error,theta_error,eq] = this.equilibriumConditions(...
                rU,q,referenceEq);
            val = [rU_error;theta_error];
            
            if any(isinf(val) | isnan(val))
                jac = nan(4,4);
                return;
            end
            
            mat = nan(4,4);
            delta = this.opt_equilibrium.TolX;
            for i = 1:4
                x_p = x;
                x_p(i) = x(i) + delta;
                rU_p = exp(x_p(1:2));
                q_p = exp(x_p(3:4));
                [rU_error_p,theta_error_p] = this.equilibriumConditions(...
                    rU_p,q_p,eq);
                mat(:,i) = [rU_error_p;theta_error_p];
            end
            jac = (mat-val)/delta; %Invalid Jacobian
            if any(isinf(jac(:)))
                jac = inf(4,4);
                val = inf(4,1);
            end
            this.previousIterEq = eq;
        end
        
        
        function eq = solveForEquilibrium(this, ...
                ini_rU,ini_q,referenceEq,findMultiplicityMode)
            verbLevel_error = 1;
            verbLevel_details = 2;
            this.previousIterEq = [];
            
            if nargin < 2 || isempty(ini_rU)
                ini_rU = [7;1.6];
            end
            if nargin < 3 || isempty(ini_q)
                ini_q = [0.55;0.6];
            end
            if nargin < 4
                referenceEq = [];
            end
            if nargin < 5
                findMultiplicityMode = false;
            end
            
            this.message(verbLevel_details, ['solveForEquilibrium: '...
                'starting with ini_rU=[' num2str(ini_rU') ...
                '], ini_q=[' num2str(ini_q') '].\n']);
            
            if this.verbosity >= verbLevel_details
                opt = optimoptions(this.opt_equilibrium, ... 
                    'display','iter');
            else
                opt = this.opt_equilibrium;
            end
            base_ini_x = [log(ini_rU);log(ini_q)];
            ini_x = base_ini_x;
            numSolutions = 0;
            attempt = 1;
            opt_rU = nan(2,0);
            opt_q = nan(2,0);
            while attempt <= this.equilibriumMaxAttempts && ...
                    (numSolutions == 0 || (findMultiplicityMode && ...
                    numSolutions == 1))
                success = false;
                if mod(attempt,2) %Odd attempts: base algorithm
                    opt = optimoptions(opt,'algorithm',this.opt_eq_base_alg);
                else
                    opt = optimoptions(opt,'algorithm',this.opt_eq_alt_alg);
                end
                try
                    [sol, fv, flag, out] = fsolve(@(x) ...
                        this.equilibriumObjFun(x,referenceEq),...
                        ini_x,opt);
                    success = max(abs(fv)) < this.max_error_equilibrium;
                catch excp                    
                    this.message(verbLevel_error, ['solveForEquilibrium: '...
                        'Error in attempt ' num2str(attempt) '/' ...
                        num2str(this.equilibriumMaxAttempts) ':\n' getReport(excp) '\n']);
                end
                if ~success                    
                    this.message(verbLevel_details, ['solveForEquilibrium: '...
                        'attempt ' num2str(attempt) ' failed.\n']);
                else
                    %Check if solution is different from previous
                    rU = exp(sol(1:2));
                    q = exp(sol(3:4));
                    differences = abs([opt_rU - rU; opt_q - q])>1e-4;
                    if numSolutions == 0 || any(differences(:))
                        numSolutions = numSolutions+1;
                        opt_rU(:,numSolutions) = rU;
                        opt_q(:,numSolutions) = q;
                        eq{numSolutions,1} = this.previousIterEq;
                        if max(abs([rU;q]-[eq{numSolutions,1}.rU;eq{numSolutions,1}.q])) > this.max_error_equilibrium
                            error('Wrong equilibrium!');
                        end
                        eq{numSolutions,1} = ...
                            this.addEquilibriumStatistics(eq{numSolutions,1});
                    end
                end
                if attempt >= 2
                    ini_x = base_ini_x + (rand(size(ini_x))-0.5)*2*2;
                end
                attempt = attempt + 1;
            end

            if numSolutions == 0
                eq = [];
                this.message(verbLevel_error, ['solveForEquilibrium: '...
                    'could not solve after ' num2str(this.equilibriumMaxAttempts) ...
                    'attempts.\n']);
                return;
            end
            if ~findMultiplicityMode
                eq = eq{1};
                this.message(verbLevel_details, ['solveForEquilibrium: '...
                    'finished with fv = ' num2str(fv') ...
                    ', opt_rU=[' num2str(opt_rU') ...
                    '], opt_q=[' num2str(opt_q') '].\n']);
            end
        end
        
        function handle = plotEquilibrium(this,eq,gGroups,vertAxisMinMax)
            if nargin < 3
                gGroups = 1:this.rhoGridSize;
            end
            greyVal = 0.6;
            styleStrings = {'-';'--';':';'-.';'--';':'};
            lineWidths = 2*[1;1;1;1;1;1];
            lineColors = [0 0 0;greyVal*ones(5,3)];
            handle = figure;
            %Formal firms
            curves(:,1) = {...
                @(z) log10(eq.fp.w_for_s(z));...
                @(z) log10(eq.fp.w_for_u(z));...
                @(z) log10(eq.fp.n_for_s(z)+eq.fp.n_for_u(z));...
                @(z) eq.fp.n_for_s(z)./(eq.fp.n_for_s(z)+eq.fp.n_for_u(z))};
            lb(1,1) = this.max2minProd;
            ub(1,1) = this.max2minProd;
            labels{1,1} = 'Formal';
            nCurves = 1;
            for g = gGroups
                if isempty(eq.fp.z_indiff{g}) && ...
                        ~eq.fp.informal_type_g(sqrt(this.max2minProd),g)
                    lb(1,1) = 1;
                    continue;
                end
                nCurves = nCurves + 1;
                curves(:,nCurves) = {...
                    @(z) log10(eq.fp.w_inf_s{g}(z));...
                    @(z) log10(eq.fp.w_inf_u{g}(z));...
                    @(z) log10(eq.fp.n_inf_s{g}(z)+eq.fp.n_inf_u{g}(z));...
                    @(z) eq.fp.n_inf_s{g}(z)./(eq.fp.n_inf_s{g}(z)+eq.fp.n_inf_u{g}(z))}; 
                lb(nCurves,1) = 1;
                if isempty(eq.fp.z_indiff{g})
                    ub(nCurves,1) = this.max2minProd;
                else
                    ub(nCurves,1) = eq.fp.z_indiff{g};
                    lb(1,1) = min(lb(1,1),ub(nCurves,1));
                end
                labels{nCurves,1} = ['Informal k=' num2str(g)];
            end
            names = {'Log_{10} wage, skilled';'Log_{10} wage, unskilled';...
                'Log_{10} employment';'Share of skilled workers'};
            for iP = 1:4
                subplot(3,2,iP,'align');
                hold on;
                for iC = 1:nCurves
                    lzv = linspace(log10(lb(iC)),log10(ub(iC)),1000)';
                    zv = 10.^lzv;
                    plot(lzv,curves{iP,iC}(zv),styleStrings{iC},...
                        'LineWidth',lineWidths(iC),...
                        'Color',lineColors(iC,:));
                end
                if iP == 4
                    legend(labels,...
                        'location','southeast');
                end
                title(names{iP});
                if nargin >= 4
                    axis([0 log10(this.max2minProd) ...
                        vertAxisMinMax(iP,1) vertAxisMinMax(iP,2)]);
                end
            end
            %Density/formal share graphs
            names = {'Firm distribution';...
                'Employment distribution'};
            labels = {'Normalized density';'Share informal'};
            styleStrings = {'--';'-'};
            lineWidths = 2*[1;1];
            lineColors = [0 0 0;greyVal*ones(5,3)];
            lzv = linspace(0,4,100000)';
            zv = 10.^lzv;
            f = this.z_density(zv);
            f_informal = zeros(size(zv));
            n_formal = zeros(size(zv));
            n_informal = zeros(size(zv));
            for g = 1:this.rhoGridSize
                f_informal = f_informal + this.rhoGridProbs(g) * ...
                    f.*eq.fp.informal_type_g(zv,g);
                n_formal = n_formal + this.rhoGridProbs(g) * ...
                    f.*(1-eq.fp.informal_type_g(zv,g)).*(...
                    eq.fp.n_for_s(zv)+eq.fp.n_for_u(zv));
                n_informal = n_informal + this.rhoGridProbs(g) * ...
                    f.*eq.fp.informal_type_g(zv,g).*(...
                    eq.fp.n_inf_s{g}(zv)+eq.fp.n_inf_u{g}(zv));
            end
            dz = diff(zv);
            meanInInterval = @(x) (x(1:(end-1))+x(2:end))/2;
            lzv_mean = log10(meanInInterval(zv));
            divideByMax = @(x) x/max(x);
            normalized = @(x) divideByMax(meanInInterval(x).*dz);
            firm_density_norm = normalized(f);
            worker_density_norm = normalized(n_formal+n_informal);
            firm_share_informal = meanInInterval(f_informal./f);
            worker_share_informal = meanInInterval(n_informal...
                ./(n_formal+n_informal));
            curves = {firm_density_norm firm_share_informal;
                worker_density_norm worker_share_informal};
            nCurves = 2;
            for iP = 1:2
                subplot(3,2,4+iP,'align');
                hold on;
                for iC = 1:nCurves                    
                    plot(lzv_mean,curves{iP,iC},styleStrings{iC},...
                        'LineWidth',lineWidths(iC),...
                        'Color',lineColors(iC,:));
                end
                if iP == 1
                    legend(labels,...
                        'location','northeast');
                end
                title(names{iP});
                xlabel('Log_{10} z');
            end
        end
       
        function newObj = clone(oldObj,numberClones)
            % produce deep copies of the model object
            if nargin==1 || numberClones==1
                props = properties(oldObj);
                newObj=Model();
                for i = 1:length(props)
                    prop=findprop(oldObj,props{i});
                    if prop.Dependent==0
                        newObj.(props{i}) = oldObj.(props{i});
                    end
                end
            else
                newObj=cell(numberClones,1);
                for i=1:numberClones
                    newObj{i}=clone(oldObj);
                end
            end
        end
      
        function message(this, verbosityLevel, str)
            paddingConstant = 3;
            if this.verbosity >= verbosityLevel
                padLeft = repmat(' ',1,paddingConstant*verbosityLevel);
                fprintf(this.verbosityOutput,[this.outputLabel ': ' padLeft str '\n']);
            end
        end
    end
end


