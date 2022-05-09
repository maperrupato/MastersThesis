%--------------------------------------------------------------------------
% Function to Set the Steady State of Markov Swithcing Model
% Marina Perrupato Mendonça
% Master's Thesis
%--------------------------------------------------------------------------
function [ss,newp,retcode] = FisLim_steadyStateFile(m, ss, params, d, id)

% Args:
%    m (rise | dsge): model object (not always needed)
%    ss (vector): endo_nbr x 1 vector of initial steady state
%    params (struct): parameter structure
%    d (struct): definitions
%    id (vector): location of the variables to calculate

% Returns:
%    - **ss** []: endo_nbr x 1 vector of updated steady state
%    - **newp** [struct]: structure containing updated parameters if any
%    - **retcode** [0|number]: return 0 if there are no problems, else return
%      any number different from 0

% - imposed(default=false): RISE computes the solution at the
%                   specified point without checking that the point solves for the
%                   steady state
% - unique (default=false): RISE computes the steady state at
%                   the ergodic distribution of the parameters. In case the
%                   probabilities are endogenous, the ergodic distribution of the
%                   parameters is itself a function of the steady state of the
%                   variables.

%% PRELIMINARY

retcode = 0;
if nargin == 1
    % list of endogenous variables to be calculated
    %----------------------------------------------
    ss = m.endogenous.name;
    
    % list of parameters herein calculated
	%-------------------------------------
	newp={};
    
else
    % no parameter to update
    %-----------------------
    newp=[];

%Economic activity
Y = 1;
C=params.C_bar;
N=params.N_bar;
A=params.A_bar;%1/N_bar;
%W_real=N_bar^(varphi)*C_bar^(siggma)/(1-tau_bar);
W_real=params.varphi_N*params.N_bar^(params.varphi)/(1-params.tau_bar); %GHH
MC=(params.epsilon-1)/(params.epsilon*params.s);

%Inflation and interest rate
Pii = 1;
Pi_D = 1;
Pi_F = 1;
Pi_x_D=1;
Pi_x_F=1;
d_D=1;
d_F=1;
R=1/params.betta;

%Fiscal Variables
G=params.G_bar; %mean IBGE data
Z=params.Z_bar; %mean IMF data
tau=(1-params.bindTax)*params.tau_bar + params.bindTax*params.tau_max; %mean IMF data
tauSdw=params.tau_bar;
delta=params.bindFL * params.delta_bar;
%T_ls=B_bar*((1/R_bar)-(1-delta)/Pi_bar)-G_bar-Z_bar+tau_bar*((C_bar+G_bar)*(1-alppha*(PSI_fbar^(eta)*Q_bar^(1-eta)))+Q_bar*Y_xbar*PSI_fbar^(eta-1));
T_ls=params.B_bar*((1/params.R_bar)-(1-delta)/params.Pi_bar)-params.G_bar-params.Z_bar+params.tau_bar*((params.C_bar+params.G_bar)*(1+params.alppha-params.alppha*params.Q_bar));
B=params.B_bar; %mean BCB data (aumentando->maior area sob dominancia fiscal)

%%Fiscal Limit
muFisLim =   params.muFisLim_Intercept;
stdFisLim  = params.stdFisLim_Intercept;
probDefault = 1/(1 + exp( params.probDefFisLim_Param_0 + params.probDefFisLim_Param_B*(B - muFisLim)));

%External Sector
M=params.alppha*(C+G);
dS=1;
TOT=1;
PSI_f=params.PSI_fbar;
Q=params.Q_bar;		
V=(params.alppha-params.alppha*params.PSI_fbar)/(1-(1/params.betta));
phi=0;
Phi_v=1;%1
Pi_x=1;
R_x=1/params.betta;
Y_x=params.Y_xbar;

%Price Setting
%x_1_D=Gamma*C^(-siggma)*Y/(1-theta_D*betta);
%x_2_D=Gamma*C^(-siggma)*Y*MC/(1-theta_D*betta);
%x_1_F=Gamma*C^(-siggma)*M/(1-theta_F*betta);
%x_2_F=Gamma*C^(-siggma)*M*PSI_f/(1-theta_F*betta);

%GHH:
x_1_D=((C+params.alpha_G*G-(params.varphi_N*N^(1+params.varphi))/(1+params.varphi))^(-params.siggma))*Y/(1-params.theta_D*params.betta);
x_2_D=((C+params.alpha_G*G-(params.varphi_N*N^(1+params.varphi))/(1+params.varphi))^(-params.siggma))*Y*MC/(1-params.theta_D*params.betta);
x_1_F=((C+params.alpha_G*G-(params.varphi_N*N^(1+params.varphi))/(1+params.varphi))^(-params.siggma))*M/(1-params.theta_F*params.betta);
x_2_F=((C+params.alpha_G*G-(params.varphi_N*N^(1+params.varphi))/(1+params.varphi))^(-params.siggma))*M*PSI_f/(1-params.theta_F*params.betta);

%Shocks
Eps_R=0;
epsilon_yx=0;
epsilon_pix=0;
epsilon_ix=0; 

%Expectational Variables
ExpA=A;
ExpG=G;    
ExpDelta=delta;	

%% Output

ys = cell2mat({arrayfun(@(x) evalin('caller', x), string(m.endogenous.name))});

% check the validity of the calculations
    %----------------------------------------
if ~utils.error.valid(ys)
    retcode=1;
else
    % push the calculations
    %----------------------
    ss(id)=ys;
end


end
