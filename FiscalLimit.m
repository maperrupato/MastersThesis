%--------------------------------------------------------------------------
% Calculating Fiscal Limit following Bi (2014,2016)
% Marina Perrupato Mendonça
% Master's Thesis
%--------------------------------------------------------------------------

%% Calculate Fiscal Limit

compute_FisLim    = true;
compute_normalfit = true;
save_FisLim       = true;
save_FisLim_fit   = true;
load_Fislim       = true;
load_Fislim_fit   = true;

reestimate_FisLim_withapproximation = false;

per = 200; 
NSimu = 150000; 
initper = 1; 

% Parametrization:
parameters_FisLim;

%%Load parameters
paramsnames = fieldnames(params);
for ii=1:length(params)
   eval([paramsnames{ii}, ' = params.', paramsnames{ii}, ';']); 
end

%Compute Fiscal Limit

if compute_FisLim
%   % enable the following two lines for parallel computing
    %parpool; 
    
    FisLim_mat = zeros(NSimu,nGrid,nGrid);
                tic

    %Iterate the state space:
    for iA = 1:nGrid %numero de pontos no grid
        for iG = 1:nGrid %numero de pontos no grid
     disp(['Iterating (A and G): ', sprintf('%d',iA), ' ', sprintf('%d',iG), ' ', sprintf('%d',1)]);
            %Simulate Fiscal Limit
            parfor iSimu = 1:NSimu
                % t = 1
                A0 = Agrid(iA);
                G0 = Ggrid(iG);
                
                Ymax = (((params.epsilon-1)/(params.epsilon*params.s))*((1-params.tau_max)/params.varphi_N))^(1/params.varphi)*((A0*params.K_bar*(G0/params.G_bar)^(params.gamma_G))^(1+(1/params.varphi)));
                Nmax = (((params.epsilon-1)/(params.epsilon*params.s))*params.K_bar*((1-params.tau_max)*A0/params.varphi_N)*(G0/params.G_bar)^(params.gamma_G))^(1/params.varphi);
                Cmax = Ymax-G0;
                Uc0 = (Cmax+(params.alpha_G*G0)-(Nmax^(1+params.varphi))/(1+params.varphi))^(-params.siggma);
                
                FisLim_mat(iSimu,iA,iG) = FisLim_mcmc(A0,G0, Uc0, per,...
                    initper, params.sigma_A, params.sigma_G, params.rho_A, params.rho_G,params.rho_Y,params.A_bar,...
                    params.G_bar,params.betta, params.epsilon, params.s, params.tau_max, params.varphi_N, params.varphi,...
                    params.siggma,params.Z_bar,params.T_LS,params.gamma_G,params.alpha_G,params.K_bar);
            end
            
        end
    end
    toc
   %parpool close
elseif load_Fislim
        load ([pathSaved filesep 'FisLim_mat.mat'],'FisLim_mat');

end

%% Normal Fit

if compute_normalfit
    
    FisLim_mu=NaN(nGrid,nGrid);
    FisLim_std=NaN(nGrid,nGrid);
    
    tic
    for iA = 1:nGrid
        for iG = 1:nGrid
            
            FisLim_vec = FisLim_mat(:,iA,iG);
            [FisLim_mu(iA,iG),FisLim_std(iA,iG)]= normfit(FisLim_vec);
            
        end
    end
    toc
elseif load_Fislim_fit
    load ([pathSaved filesep 'FisLim_normfit.mat'],'FisLim_mu','FisLim_std');
end


%% Save Results
pathSaved = ['C:/Users/marin/OneDrive/Documents/Mestrado/PaperVerao/Saved'];
%PC Puc:
%pathSaved = ['C:/Users/marina.mendonca/OneDrive/Documents/Mestrado/PaperVerao/Saved'];

if save_FisLim
    save([pathSaved filesep 'FisLim_mat.mat'],'FisLim_mat');
end

if save_FisLim_fit
    save([pathSaved filesep 'FisLim_normfit.mat'],'FisLim_mu','FisLim_std');
end


%% Graphs

pdfDist_L = makedist('Normal', FisLim_mu(1,1) /(4*params.Y_bar), FisLim_std(1,1)/(4*params.Y_bar));
pdfDist_M = makedist('Normal', FisLim_mu(5,5) /(4*params.Y_bar), FisLim_std(5,5)/(4*params.Y_bar));
pdfDist_H = makedist('Normal', FisLim_mu(11,11) /(4*params.Y_bar), FisLim_std(11,11)/(4*params.Y_bar));
pdfDist_better = makedist('Normal', FisLim_mu(1,11) /(4*params.Y_bar), FisLim_std(1,11)/(4*params.Y_bar));
pdfDist_SS = makedist('Normal', FisLim_mu(6,6) /(4*params.Y_bar), FisLim_std(6,6)/(4*params.Y_bar));

x = 0.0:0.001:1.75;
figure('Name','Fiscal Limit CDF - Marina Calibration')
color_FL=[0.27 0.584 0.576
    0.76	0.37 0.75
    0.89 0.34 0.68
    0.15	0.55	0.88
    0.74	0.48	0.65];
colororder(color_FL)
plot(x*100, cdf(pdfDist_L, x)*100,'linewidth',2);
hold on
plot(x*100, cdf(pdfDist_M, x)*100,'linewidth',2);
hold on
plot(x*100, cdf(pdfDist_H, x)*100,'linewidth',2);
hold on
plot(x*100, cdf(pdfDist_better, x)*100,'linewidth',2);
hold on
plot(x*100, cdf(pdfDist_SS, x)*100,'linewidth',2);
xticks([0:25:175])
yticks([0:10:100])
xlabel('Debt-GDP (%)','fontsize',10)
legend({'(A,G) Low', '(A,G) Medium','(A,G) High', 'Lowest \mu', 'Steady State \mu'},'Location','northwest')
title('Fiscal Limit CDF')

x = 0.0:0.001:1.5;
figure('Name','Fiscal Limit CDF - Steady State')
color_FLSS=[0.74	0.48	0.65];
colororder(color_FLSS)
plot(x*100, cdf(pdfDist_SS, x)*100,'linewidth',2);
xticks([0:25:150])
xlabel('Debt-GDP (%)','fontsize',10)
legend({'Steady State \mu'},'Location','northwest')
title('Fiscal Limit CDF')

%% Fiscal Limit Approximations and their parameters

% Approximate Fiscal Limit's mu and std by a Polynomial Function
FisLim_polyapprox;