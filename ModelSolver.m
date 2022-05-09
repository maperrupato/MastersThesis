%--------------------------------------------------------------------------
% Endogenous Markov Switching Model Solver
% Marina Perrupato Mendonça
% Master's Thesis
%--------------------------------------------------------------------------
%% load the RISE paths
addpath('C:\Users\marin\OneDrive\Documents\Mestrado\PaperVerao\RISE_toolbox-master')

% Load RISE. Run this only once.
rise_startup();

%% RISE the model
m = rise('ModelOpenBi_subsidy', 'steady_state_file', 'FisLim_steadyStateFile','parameters',params);

regimeNames     = { 'Tax rate below max / Debt below fiscal limit', ...
    'Tax rate at max / Debt below fiscal limit', ...
    'Tax rate below max / Debt crossed fiscal limit', ...
    'Tax rate at max / Debt crossed fiscal limit'};

%Parameterize, solve the model, and assign the data (Policy Rules)
solveOrder  = 1;
solve_derivatives_type = 'symbolic'; % symbolic, automatic, numerical  %%% in case of optimal rule, cannot be symbolic
debugMode   = false; % true, false
steady_state_imposed = true;  % if imposed, then does not check whether it is a steady state (solves around the reference regime ss)
steady_state_unique = false; % the steady state is unique: solves around the "ergodic mean" and does not run the steady state file

m=solve(m,'debug', debugMode, ...
    'solve_order', solveOrder,...    % may not work under endogenous regime switching
    'steady_state_imposed', steady_state_imposed, ...   % if imposed, then does not check whether it is a steady state
    'steady_state_unique',  steady_state_unique);      % the steady state is unique: solves around the ergodic mean and does not run the steady state file

cfm=is_stable_system(m, 'stability_algorithm', 'cfm');
gmh=is_stable_system(m, 'stability_algorithm', 'gmh');

%% Steady state

m.markov_chains.regimes

ssMatrix     = [];
ssVarNames   = [];

nRegimes     = m.markov_chains.regimes_number;
ssMatrix     = [ssMatrix, [m.solution.ss{:}]];
ssVarNames   = [ssVarNames; strcat("m", "; Regime ", string(1:nRegimes)')];

tSS = array2table(ssMatrix);
tSS.Properties.VariableNames    = ssVarNames;
tSS.Properties.RowNames         = m.endogenous.name;
disp(tSS);



%% IRF graphs
%creating vectors of variables and shocks
varlist={'Y','C','M','Pi_D','Pi_F','Pii','dS', 'Q', 'TOT','R','G','N','V','B','tau'};
varnames={'Output','Consumption','Imports','Inflation - Domestic Prod.','Inflation - Imported',...
    'Inflation','\Delta Nominal Exchange Rate','Real Exchange Rate','Terms of Trade','Gross Interest Rate',...
    'Government Expenses','Labor','Net foreign asset position','Government Debt','Income Tax'};
shocklist={'e_G','e_R'};
shocknames={'Fiscal Shock','Monetary Shock'};
shockSize= 2;
simul_honor_constraints                 = true;
simul_honor_constraints_through_switch  = true;

myirf=irf(m,'irf_shock_list',shocklist,'irf_var_list',varlist,'irf_periods',40,...
    'simul_honor_constraints', simul_honor_constraints, ...
    'irf_shock_sign', shockSize,...
    'simul_honor_constraints_through_switch',simul_honor_constraints_through_switch);

for ishock=1:numel(shocklist)
    shock=shocklist{ishock};
    shockname=shocknames{ishock};
    figure('name',['IRF to ', shock])
    color_irf=[0.00 0.50 0.50
        0.27	0.51 0.66
        0.004 0.70 0.69
        0.16 0.29 0.37];
    for ivar=1:numel(varlist)
        colororder(color_irf)
        subplot(5,3,ivar)
        plot(myirf.(shock).(varlist{ivar}),'linewidth',2)
        %xticks([0:10:40]);
        title(varnames{ivar})
    end
    % add a bit space to the figure
    fig = gcf;
    %fig.Position(3) = fig.Position(3) + 100;%250;
    % add legend
    Lgnd = legend({'Regime 1','Regime 2','Regime 3','Regime 4'},'Orientation','horizontal');
    Lgnd.Position(1) = 0.18;%0.9;
    Lgnd.Position(2) = 0.02;%0.5;
    [~,tmp]=sup_label([shockname],'t');
    set(tmp,'fontsize',15)
end

regimeNames  = {'Tax rate below max / Debt below fiscal limit', ...
                   'Tax rate at max / Debt below fiscal limit', ...
                   'Tax rate below max / Debt crossed fiscal limit', ...
                   'Tax rate at max / Debt crossed fiscal limit'};


%Regime 1
for ishock=1:numel(shocklist)
    shock=shocklist{ishock};
    shockname=shocknames{ishock};
    figure('name',['IRF to ', shock])
    for ivar=1:numel(varlist)
        subplot(5,3,ivar)
        plot(myirf.(shock).(varlist{ivar}).data(:,1),'linewidth',2,'Color',[0.00 0.50 0.50])
        title(varnames{ivar})
    end
    [~,tmp]=sup_label([shockname],'t');
    set(tmp,'fontsize',15)
end

%Regime 2
for ishock=1:numel(shocklist)
    shock=shocklist{ishock};
    shockname=shocknames{ishock};
    figure('name',['IRF to ', shock])
    for ivar=1:numel(varlist)
        subplot(5,3,ivar)
        plot(myirf.(shock).(varlist{ivar}).data(:,2),'linewidth',2,'Color',[0.27 0.51 0.66])
        title(varnames{ivar})
    end
    [~,tmp]=sup_label([shockname],'t');
    set(tmp,'fontsize',15)
end