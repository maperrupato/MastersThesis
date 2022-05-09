%--------------------------------------------------------------------------
% Approximate Fiscal Limit's mu and std by a Polynomial Function
% Marina Perrupato Mendonça
% Master's Thesis
%--------------------------------------------------------------------------

%% Approximate Fiscal Limit by a Polynomial Function

if compute_normalfit || load_Fislim_fit
    
    % Tauchen Method - A and G
    %Common information:
    nGrid=11;
    nStd=3;
    distName='Normal';
    valMean=0;
    
    %Discretize the technology (A) shock:
    [agrid,aprob]=ARtoMarkov(nGrid,nStd,distName,valMean,params.sigma_A,params.rho_A);
    %Obtaining the grid for A in level:
    Agrid=params.A_bar*exp(agrid);
    
    %Discretize the government spending (G) shock:
    [ggrid,gprob]=ARtoMarkov(nGrid,nStd,distName,valMean,params.sigma_G,params.rho_G);
    %Obtaining the grid for G in level:
    Ggrid=params.G_bar*exp(ggrid);
else
    nGrid = 1;
    Agrid = 0;
    Ggrid = 0;
end

% Build grid with the discretizations' node
nodeGrid  = [(1:length(Agrid))',(1:length(Agrid))'];
nGridComb = allcomb(nodeGrid(:,1),nodeGrid(:,2));

% Build grid to generate mu values according to shocks
%muGrid = [Agrid, Ggrid];
muGrid = [agrid, ggrid];
muComb = allcomb(muGrid(:,1), muGrid(:,2));

% Build grid to generate std values according to shocks
%stdGrid = [Agrid, Ggrid];
stdGrid = [agrid, ggrid];
stdComb = allcomb(stdGrid(:,1), stdGrid(:,2));

% Create functions to return mu and std according to A and G
fRead_mu = @(x) FisLim_mu(nGridComb(x,1),nGridComb(x,2));
fVals_mu = arrayfun(fRead_mu,(1:size(nGridComb,1))');

fRead_std = @(x) FisLim_std(nGridComb(x,1),nGridComb(x,2));
fVals_std = arrayfun(fRead_std, (1:size(nGridComb,1))');

% Approximate a polynomial function for the fiscal limit mu and std
maxPower = 1;
global fitModel_mu
global fitModel_std

if nGrid == 1
    % Fixed point function
    fitModel_mu     = @(x1,x2) fVals_mu;
    fitModel_std    = @(x1,x2) fVals_std;
elseif isempty(muComb)
    % Fixed point function
    fitModel_mu     = @(x1,x2) fVals_mu;
    fitModel_std    = @(x1,x2) fVals_std;
else
    % Obtaining the matrix whithout repetition
    [muComb, idxUn] = unique(muComb, 'rows');
    fVals_mu    = fVals_mu(idxUn);
    
    [stdComb, idxUn] = unique(stdComb, 'rows');
    fVals_std    = fVals_std(idxUn);
    
    % Multivariate Polinomyal Regression - mu
    sRegMu = MultiPolyRegress(muComb, fVals_mu, maxPower);
    fitModel_mu     = @(x1,x2) sRegMu.PolynomialExpression(x1,x2);
    
    %fMu_str = ['@(x1,x2) ' num2str(sRegMu.Coefficients(1)) ...
    %           '+' num2str(sRegMu.Coefficients(2)) '*x1 '...
    %      '+' num2str(sRegMu.Coefficients(3)) '*x2 '...
    %      ];
    %fitModel_mu = str2func(fMu_str);
    
    % Multivariate Polinomyal Regression - std
    sRegStd = MultiPolyRegress(stdComb, fVals_std, maxPower);
    fitModel_std     = @(x1,x2) sRegStd.PolynomialExpression(x1,x2);
    
    %fStd_str = ['@(x1,x2) ' num2str(sRegStd.Coefficients(1)) ...
    %           '+' num2str(sRegStd.Coefficients(2)) '*x1 '...
    %      '+' num2str(sRegStd.Coefficients(3)) '*x2 '...
    %      ];
    %fitModel_std = str2func(fStd_str);
end

% Approximation of the default probability to a logistic function
if compute_normalfit == true
    if reestimate_FisLim_withapproximation
        % Approximate with a logistic function
         FisLim_logisticapprox;
    else
        load([pathSaved filesep 'FisLim_Logisticapprox.mat']);
    end
end

%% Parametrize Probabilities Functions

if sum(muGrid(:)) + sum(stdGrid(:)) == 0 %|| compute_normalfit
    
    % Fiscal Limit Distribution - Normal Function
    %mu
    params.muFisLim_Intercept = fVals_mu;
    params.muFisLim_A         = 0;
    params.muFisLim_G         = 0;
    %std
    params.stdFisLim_Intercept = fVals_std;
    params.stdFisLim_A         = 0;
    params.stdFisLim_G         = 0;
    
    % Default Probability - Logistic Function
    params.probDefFisLim_Param_0  = 0;
    params.probDefFisLim_Param_B  = 0;
    params.probDefFisLim_Param_A  = 0;
    params.probDefFisLim_Param_G  = 0;
else
    
    % Fiscal Limit Distribution - Normal Function
    %mu
    params.muFisLim_Intercept = sRegMu.Coefficients(1);
    params.muFisLim_A         = sRegMu.Coefficients(2);
    params.muFisLim_G         = sRegMu.Coefficients(3);
    %std
    params.stdFisLim_Intercept = sRegStd.Coefficients(1);
    params.stdFisLim_A         = sRegStd.Coefficients(2);
    params.stdFisLim_G         = sRegStd.Coefficients(3);
    
    % Default Probability - Logistic Function
    params.probDefFisLim_Param_0 = probFisLim_logistic.Coefficients.Estimate(1);
    params.probDefFisLim_Param_B = probFisLim_logistic.Coefficients.Estimate(2);
    params.probDefFisLim_Param_A = probFisLim_logistic.Coefficients.Estimate(3);
    params.probDefFisLim_Param_G = probFisLim_logistic.Coefficients.Estimate(4);
end


















