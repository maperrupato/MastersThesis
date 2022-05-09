%--------------------------------------------------------------------------
% Fit Default Probability to a Logistic Function
% Marina Perrupato Mendonça
% Master's Thesis
%--------------------------------------------------------------------------

% Set-up simulation
nFisLimSims       = 1000;
%rng(11) %coef_G>0 - erro pequeno 5e-05
%rng(9)  %coef_G<0 - 3e-04
rng(68)  %coef_G<0 - 1.2e-04
randList          = randn(nFisLimSims);

% Get data
stateComb  = muComb;
muFisLim   = fVals_mu;
stdFisLim  = fVals_std;
nFisLim    = length(muFisLim);

% Build Fiscal Limit sample
sampleFisLim = [];
for iComb = 1:nFisLim
    distFisLim    = muFisLim(iComb) + stdFisLim(iComb)*randn(nFisLimSims, 1);
    for iSorted = 1:nFisLimSims
        probFisLim = normcdf(distFisLim, muFisLim(iComb), stdFisLim(iComb));
    end
    
    %sampleFisLim(:,1): probabilities
    %sampleFisLim(:,2): mu found in distFisLim
    %sampleFisLim(:,3): mu found in fvals_mu - mu corresp to A and G in stateComb
    %sampleFisLim(:,4 and 5): A and G combinations
    
    sampleFisLim = [sampleFisLim; [probFisLim, distFisLim, repmat(muFisLim(iComb), nFisLimSims, 1),...
    repmat(stateComb(iComb,:), nFisLimSims, 1)]];
end

% Fit non-linear model
y = sampleFisLim(:,1); 
X = sampleFisLim(:,2:end);    
logistfun = @(b,x) 1./(1 + exp( b(1) + b(2)*(x(:,1) - x(:,2)) + b(3)*x(:,3) + b(4)*x(:,4)));
beta0 = [0.5 0.5 0.5 0.5];
%beta0 = [0 -10 0 0];

probFisLim_logistic = fitnlm(X,y,logistfun,beta0);
%mdlFL_logistic.disp

% Compare values
comparisonValues = [logistfun(probFisLim_logistic.Coefficients.Estimate,X), fFiscalLimit_Prob(X(:,1), X(:,3),X(:,4))];
diff = comparisonValues(:,1)-comparisonValues(:,2);
disp(['Mean difference of the fiscal limit approximation with a logistic function:' num2str(mean(diff))]);

%% Save
save([pathSaved filesep 'FisLim_Logisticapprox.mat'], 'probFisLim_logistic');

