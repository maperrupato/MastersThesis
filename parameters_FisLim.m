%--------------------------------------------------------------------------
% Parameters for Fiscal Limit Calculation
% Marina Perrupato Mendonça
% Master's Thesis
%--------------------------------------------------------------------------

%% Parametrization:
%Preferences
params.alppha=0.25;
params.betta=0.987; %0.989, 0.987 norm: 0.982
params.varphi=1; %Jonas:1.6, Eduardo: 1
params.siggma=2.132; %Eduado: 2.132 (Simstruct.params), 1.3 (script)/ Carlos: 1.3
params.eta=1;%Jonas:0.4501
params.chi=0.03;
params.epsilon=5; %Eduardo: 11 Marina: 5
params.alpha_G=0;%0.54204 (Simstruct.params) - aumenta espaço fiscal
params.gamma_G=0;%0.16301 / 0.3 - reduz espaço fiscal

%Price indexation
params.delta_D=0.65; %Jonas:0.5511
params.delta_F=0.65; %Jonas:0.935

%Price rigidity
params.theta_D=0.65; %Jonas:0.8713
params.theta_F=0.65; %Jonas:0.1339

%Fiscal policy
params.rho_tau=0.862; %Eduardo:0.862, OLS:0.85553, Bi(2017): 0.8 (culpado por reduzir a área dominancia fiscal)
params.rho_B=0.108; %Bi (2017)=0.047,before: estimated by OLS=0.1083
params.delta_bar=0.05; %Fiscal limit never bind-> delta=0, for all t
params.s=(params.epsilon-1)/params.epsilon; %Marina (Preston): (params.epsilon-1)/params.epsilon   Eduardo:1

%Regime Switching
%%Fiscal Limit is reached or not
params.bindFL_FL_1 = 0;
params.bindFL_FL_2 = 1;

%%Peak of Laffer Curve is reached or not
params.bindTax_Tax_1 = 0;
params.bindTax_Tax_2 = 1;

%Monetary policy
params.lambda_pi=1.99; %Jonas: 1.52, Carlos 1.99
params.lambda_y=0.2;%0.5; %quando 0 recupera Leeper (1991)/%Jonas:0.1279
params.lambda_s=0.31; %0.31, Jonas:4.9353,menor é 2.32 %Carlos: 0.31
params.sigma_R=0.79; %Jonas:0.335

%Steady State
params.Y_bar=1;
params.Pi_bar=1;
params.R_bar=1/params.betta;
params.G_bar=0.216; %0.186; %mean IBGE data: 0.186, Eduardo: 0.206
params.Z_bar=0.142; %mean IMF data
params.tau_bar=0.391; %mean IMF data: 0.3937 / Eduardo: 0.391
params.B_bar=2.58; %mean BCB data: 2.58 (aumentando->maior area sob dominancia fiscal)/ Eduardo: 2.476
params.tau_max=(params.varphi/(1+params.varphi));

params.PSI_fbar=(params.epsilon-1)/(params.epsilon*params.s);
params.Q_bar=(params.epsilon-1)/(params.epsilon*params.s);

params.C_bar=params.Y_bar-params.G_bar;
%params.N_bar=(((params.epsilon-1)/params.epsilon)*(1/params.s)*(1-tau_bar))^(1/(1+params.varphi));
params.N_bar=1/3;
params.K_bar=18;
params.A_bar=1/(params.N_bar*params.K_bar);

params.varphi_N=(1-params.tau_bar)/(params.N_bar^(1+params.varphi));
params.Y_xbar=params.alppha*(params.Q_bar^(-params.eta))*(params.C_bar+params.G_bar);
params.Pi_xbar=1;
params.R_xbar=1/params.betta;

params.T_LS=params.B_bar*((1/params.R_bar)-(1-params.delta_bar)/params.Pi_bar)-params.G_bar-params.Z_bar+params.tau_bar*((params.C_bar+params.G_bar)*(1+params.alppha-params.alppha*params.Q_bar));
%params.T_LS=0;

%Shocks
params.rho_A=0.9332; %Jonas:0.3485/ Carlos:0.5/ Eduardo:0.9332(Simstruct.params) 0.964 (script)
params.rho_G=0.7949; %Bi (2017),before: Eduardo=0.7949(Simstruct.params) 0.951 (script)
params.rho_Y=0.132; %Eduardo: atual 0.13201 (Simstruct.params) 0.036(script)	
params.rho_phi=0.50;%Preston:0.8/%Jonas: 0.629
params.sigma_A=0.0046; %Eduardo:0.0046(Simstruct.params),0.007 (script) / Carlos:1    
params.sigma_phi=1;%Preston:0.5
params.sigma_G=0.013; %Eduardo=0.013(Simstruct.params), 0.012 (script)

%VAR Foreign Economy
%Calibração Jonas
params.a0_py = -0.5251; 
params.a0_iy = -0.2648; 
params.a0_ip = 0.0103;
params.a1_yy = 0.9038;
params.a1_yp = 0.0839;
params.a1_yi = -0.0114;
params.a1_py = -0.3588;
params.a1_pp = -0.1458;
params.a1_pi = 0.0444;
params.a1_iy = -0.1765;
params.a1_ip = -0.0155;
params.a1_ii = -0.9683;
params.sigma_yx = 0.0034;
params.sigma_pix = 0.0062;
params.sigma_ix = 0.0048;

%% Tauchen Method (1986)

%Common information:
nGrid=11;
nStd=3;
distName='Normal';
valMean=0;

%Discretize the technology (A) shock:
%Equation: log(A/A_bar)=rho_A*log(A(-1)/A_bar)+sigma_A*e_A
[agrid,aprob]=ARtoMarkov(nGrid,nStd,distName,valMean,params.sigma_A,params.rho_A);
%Obtaining the grid for A in level:
Agrid=params.A_bar*exp(agrid);

%Discretize the government spending (G) shock:
%Equation: log(G/G_bar)=rho_G*log(G(-1)/G_bar)+rho_Y*log(Y(-1)/Y_bar)+sigma_G*e_G;
[ggrid,gprob]=ARtoMarkov(nGrid,nStd,distName,valMean,params.sigma_G,params.rho_G);
%Obtaining the grid for G in level:
Ggrid=params.G_bar*exp(ggrid);



