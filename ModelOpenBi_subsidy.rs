%--------------------------------------------------------------------------
% DSGE Open Economy with Fiscal Limits - CRRA or GHH Utility Function
% Marina Perrupato Mendonça
% Master's Thesis
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Declare Endogenous Variables
%--------------------------------------------------------------------------
endogenous

C	    "$C$"        %'Consumption',				
W_real	"$W^{real}$" %'Real Wage',
N	    "$N$"        %'Hours worked',	
MC	    "$MC$"       %'Real marginal costs',				
Y	    "$Y$"        %'Output',				 
R	    "$R$"        %'Nominal Interest Rate',				
Pii	    "$\Pii$"     %'Inflation',
V       "$V$"        %'Net foreign asset position',
Phi_v   "$\Phi(V)$"  %'Risk premium function',
M       "$M$"        %'Import',

Pi_x_D  "$\Pi_{D}^{x}$" %'Readjusted domestic firms inflation',
x_1_D   "$x_{1,D}$"     %'auxiliar variable domestic firms price setting',
x_2_D   "$x_{2,D}$"     %'auxiliar variable domestic firms price setting',
Pi_D	"$\Pi_{D}$"     %'Domestic goods Inflation',
d_D     "$d_D$"         %'Domestic prices dispersion index',

Pi_F	"$\Pi_{F}$"      %'Imported goods Inflation',
Pi_x_F  "$\Pi_{F}^{x}$"  %'Readjusted imported firms inflation',
x_1_F   "$x_{1,F}$"      %'auxiliar variable importerd firms price setting',
x_2_F   "$x_{2,F}$"      %'auxiliar variable imported firms price setting',
d_F     "$d_{F}$"        %'Import prices dispersion index',

%External sector
Q	   "$Q$"        %'Real exchange rate',
dS     "$\Delta S$" %'Variation of Nominal exchange rate',
TOT	   "$ToT$"      %'Terms of Trade',
PSI_f  "$\Psi_f$"   %'Deviation from Law of One Price',
		
%%Foreign Economy
Pi_x    "$\Pi^{*}$"  %'Foreign Inflation',
Y_x     "$Y^{*}$"    %'Foreign output',
R_x     "$R^{*}$"    %'Foreign interest rate',
						
%	FISCAL	VARIABLES			
B	   "$B$"         %'Real Debt',				
tau	   "$\tau$"      %'Income tax rate',
tauSdw "$tau^{sdw}$" %'Income tax shadow rate',				
G	   "$G$"         %'Government expenses',
Z      "$Z$"         %'Lump sum transfers',
T_ls   "$T^{LS}$"    %'Steady state lump-sum tax',
delta  "$\delta$"    %'Default rate',

muFisLim    "$E_t \mathcal{B}_t$" %'Expected Fiscal Limit mu ~N(muFisLim,stdFisLim)',
stdFisLim   "$\sigma_t(\mathcal{B}_t)$" %'Expected Standard Dev. Fiscal Limit~N(muFisLim,stdFisLim)',
probDefault "$P_t(B_t > \mathcal{B}^{*}_{t+1})$" %'Expected Default Probability - prob debt exceed FisLim',
%polDef      "$\mathcal{D}_t$" %(COLOCOU NA TFP, ATÉ O MOMENTO NÃO FAREI ESSA MUDANÇA)

%	SHOCKS	AR(1)	VARIABLES		
A	   "$A$"               %'AR(1) technology process',				
Eps_R  "$\varepsilon^{R}$" %'AR(1) monetary policy shock process',
phi    "$\phi$"            %'AR(1) risk premium shock process',
epsilon_yx  "$\varepsilon_{y}^{*}$"   %'external demand shock',
epsilon_pix "$\varepsilon_{\pi}^{*}$" %'international inflation shock',
epsilon_ix  "$\varepsilon_{i}^{*}$"   %'international monetary shock'

% EXPECTATION VARIABLES
ExpA     "$E_t A_{t+1}$"      %'Expected technology process',	
ExpG     "$E_t G_{t+1}$"      %'Expected government expenses',	
ExpDelta "$E_t \delta_{t+1}$" %'Expected ',	

%--------------------------------------------------------------------------
% Define shock variables
%--------------------------------------------------------------------------
exogenous

e_A	        "$\epsilon_{A}$"          %'Technology shock',				
e_R	        "$\epsilon_{R}$"          %'Monetary policy shock',
e_phi       "$\epsilon_{\phi}$"       %'Risk premium shock',
e_G         "$\epsilon_{G}$"          %'Fiscal policy shock',
e_yx        "$\epsilon_{y}^{*}$"      %'external demand shock',
e_pix       "$\epsilon_{\pi}^{*}$"    %'international inflation shock',
e_ix        "$\epsilon_{i}^{*}$"      %'international monetary shock'

%--------------------------------------------------------------------------					
%	Define	Parameters			
%--------------------------------------------------------------------------				
parameters 
betta	  "$\beta$"      %'Discount factor',	
siggma	  "$\sigma$"     %'Inverse EIS',			
varphi	  "$\varphi$"    %'Inverse Frisch elasticity',
varphi_N "$\varphi_N$"  %'Disutility of Labor',
chi       "$\chi$"       %'Elasticity of risk premium to external indebtedness',
epsilon	  "$\epsilon$"   %'Demand elasticity',			
delta_D   "$\delta_{D}$" %'Domestic indexation',
theta_D   "$\theta_{D}$" %'Domestic price rigidity',
delta_F   "$\delta_{F}$" %'Imported indexation',
theta_F   "$\theta_{F}$" %'Imported price rigidity',
eta	      "$\eta$"       %'Substitution elasticity between domestic and imported',				
alppha    "$\alpha$"     %'Economic openmess',		
alpha_G   "$\alpha_G$"   %'Subs gov and private consumption',
gamma_G   "$\gamma_G$"   %'Capital enhancement effect',

%Monetary Policy
lambda_y  "$\lambda_{y}$"   %'Output feedback Taylor Rule',	
lambda_pi "$\lambda_{\pi}$" %'Inflation feedback Taylor Rule',
lambda_s  "$\lambda_{s}$"   %'Devaluation feedback Taylor Rule',

%Fiscal Policy
rho_tau   "$\rho_{\tau}$"   %'Autocorrelation Income tax',
rho_B     "$\rho_{B}$"      %'Income tax response to debt',
rho_G	  "$\rho_{G}$"      %'Autocorrelation fiscal shock',			
rho_Y     "$\rho_{Y}$"      %'Expense response to product',
sigma_G   "$\sigma_{G}$"    %'Government expenses shock variance',
s         "$s$"            %'Subsidy rate to labor and imports',
delta_bar "$\bar{\delta}$"  %'Default rate',
tau_max    "$\tau^{max}$"    %'Income Tax at Peak of Laffer Curve',

%%Fiscal Limit Distribution - Normal Function
% Multivariate Polinomyal Regression- mu
muFisLim_Intercept "$\beta^{\mu}_{0}$"
muFisLim_A         "$\beta^{\mu}_{A}$"
muFisLim_G         "$\beta^{\mu}_{G}$"

% Multivariate Polinomyal Regression- std
stdFisLim_Intercept "$\beta^{\std}_{0}$"
stdFisLim_A         "$\beta^{\std}_{A}$"
stdFisLim_G         "$\beta^{\std}_{G}$"

%%Default Probability - Logistic Function Approximation
probDefFisLim_Param_0  "$\gamma_{0}$"
probDefFisLim_Param_B  "$\gamma_{B}$"
probDefFisLim_Param_A  "$\gamma_{A}$"
probDefFisLim_Param_G  "$\gamma_{G}$"

%Shocks
rho_A	    "$\rho_{A}$"           %'Autocorrelation technology shock',								
rho_phi	    "$\rho_{\phi}$"        %'Autocorrelation risk premium shock',
sigma_A	    "$\sigma_{A}$"         %'Technology shock variance',				
sigma_R	    "$\sigma_{R}$"         %'Monetary policy shock variance',				
sigma_phi	"$\sigma_{\phi}$"      %'Risk premium shock variance',	

%Steady State parameters
Y_bar      "$\bar{Y}$"         %'Steady state output',
C_bar      "$\bar{C}$"         %'Steady state consumption',
N_bar      "$\bar{N}$"	       %'Steady state labor',
K_bar      "$\bar{K}$"         %'Capital Stock',
A_bar      "$\bar{A}$"         %'Steady state TFP',	
Pi_bar     "$\bar{\Pi}$"       %'Steady state inflation',	
G_bar      "$\bar{G}$"         %'Steady state government expenses',	
tau_bar    "$\bar{\tau}$"      %'Steady state income tax',
B_bar      "$\bar{B}$"         %'Steady state real debt',
Z_bar      "$\bar{Z}$"         %'Steady state lump sum transfers',
R_bar      "$\bar{R}$"         %'Steady state interest rate',
Pi_xbar    "$\bar{\Pi^{*}}$"   %'Steady state external inflation',
R_xbar     "$\bar{R^{*}}$"     %'Steady state external interest rate',
Y_xbar     "$\bar{R^{*}}$"     %'Steady state external output'
PSI_fbar   "$\bar{\Psi}_{f}$"  %'Steady state LOP deviation',
Q_bar      "$\bar{Q}$"         %'Steady state real exchange rate',
T_LS

%Foreign Economy VAR
a0_py     "$a_{0,\pi y}$"
a0_iy     "$a_{0,i y}$"
a0_ip     "$a_{0,i\pi}$"
a1_yy     "$a_{1,yy}$"
a1_yp     "$a_{1,y\pi}$"
a1_yi     "$a_{1,y i}$"
a1_py     "$a_{1,\pi y}$"
a1_pp     "$a_{1,\pi \pi}$"
a1_pi     "$a_{1,\pi i}$"
a1_iy     "$a_{1,i y}$"
a1_ip     "$a_{1,i \pi}$"
a1_ii     "$a_{1,i i}$"
sigma_yx  "$\sigma_{y}^{*}$"
sigma_pix "$\sigma_{\pi}^{*}$"
sigma_ix  "$\sigma_{i}^{*}$"

% SWITCHING PARAMETERS SHOULD NOT BE DECLARED HERE!
%bindFL  "$\mathbbm{1}_{\text{fiscal limit}}$"  %'If fiscal limit is bind =1',
%bindTax "$\mathbbm{1}_{\tau}^{max}$"           %'If income tax is at max =1',

%--------------------------------------------------------------------------					
% Markov Switching Structure (4 regimes) 
%--------------------------------------------------------------------------								

%Markov Chain FL with 2 states - Fiscal Limit is reached or not
parameters(FL,2) 
 bindFL  "$\mathbbm{1}_{\text{fiscal limit}}$" %'If fiscal limit is bind =1',
%parameters FL_tp_1_2 FL_tp_2_1 %Exogenous probabilities

%Markov Chain Tax with 2 states - Peak of Laffer Curve is reached or not
parameters(Tax,2)
 bindTax "$\mathbbm{1}_{\tau}^{max}$"       %'If income tax is at max =1',

%--------------------------------------------------------------------------					
%	Model				
%--------------------------------------------------------------------------								
model

% Endogenous Switching Probabilities
 %Default Probability - Current Debt Exceeds Fiscal Limit
   ! FL_tp_1_2 = 1/(1 + exp( probDefFisLim_Param_0 + probDefFisLim_Param_B*(B - muFisLim) +...
                 probDefFisLim_Param_A*(ExpA-A_bar) + probDefFisLim_Param_G*(ExpG-G_bar))) ;
   ! FL_tp_2_1 = 1;

 %Peak of the Laffer Curve
   ! Tax_tp_1_2 = (tauSdw > tau_max) ; %0.5 is the maximum tax rate
   ! Tax_tp_2_1 = (tauSdw < tau_max) ;

%[name='Expected fiscal limit mean']
muFisLim  = muFisLim_Intercept + muFisLim_A*(ExpA-A_bar) + muFisLim_G*(ExpG-G_bar) ;

%[name='Expected fiscal limit standard deviation']
stdFisLim = stdFisLim_Intercept + stdFisLim_A*(ExpA-A_bar) + stdFisLim_G*(ExpG-G_bar) ;

%[name='Expected default probability']
probDefault = 1/(1 + exp( probDefFisLim_Param_0 + probDefFisLim_Param_B*(B - muFisLim) +...
              probDefFisLim_Param_A*(ExpA-A_bar) + probDefFisLim_Param_G*(ExpG-G_bar) )) ;

%[name='1. Labor Supply']		
%W_real=(N^(varphi)*C^(siggma))/(1-tau);
W_real=(varphi_N*N^(varphi))/(1-tau); %GHH

%[name='2. Euler Equation']		
%1/R=betta*(1-delta)*((Gama(+1)/Gama)*(1/Pii(+1))*(C/C(+1))^(siggma));
1/R=betta*(1-delta)*((1/Pii(+1))*((C+alpha_G*G-(varphi_N*N^(1+varphi))/(1+varphi))/(C(+1)+...
   alpha_G*G(+1)-(varphi_N*N(+1)^(1+varphi))/(1+varphi)))^(siggma)); %GHH

%[name='3. Uncovered Interest Parity']	
%Gama(+1)*(1/Pii(+1))*(C(+1)^(-siggma))*(R_x*Phi_v*dS(+1)-(1-delta)*R)=0;
(1/Pii(+1))*(C(+1)+alpha_G*G(+1)-(varphi_N*N(+1)^(1+varphi))/(1+varphi))^(siggma)*(R_x*Phi_v*dS(+1)-(1-delta(+1))*R)=0;

%[name='4. Risk premium']
Phi_v=exp(-chi*V+phi);

%[name='5. Monetary policy']
R=(1/betta)*((Pii/Pi_bar)^(lambda_pi))*((Y/Y_bar)^(lambda_y))*dS^(lambda_s)*exp( Eps_R);

%[name='6. Real debt dynamics']
B=R*((1-delta)*(B(-1)/Pii)+G+Z+T_ls-tau*(((Q*Y)/(PSI_f*TOT))+(C+G)*(1-...
  ((Q/(PSI_f*TOT))^(1-eta))*(1-alppha+(alppha*PSI_f*TOT^(1-eta))))));

%[name='7. Income tax']
%log(tau/tau_bar)=rho_tau*log(tau(-1)/tau_bar)+rho_B*log(B(-1)/B_bar));
%log(tau)=(1-bindTax)*(log(tau_bar)+rho_tau*log(tau(-1)/tau_bar)+rho_B*log(B(-1)/B_bar))+bindTax*(log(tau_max));
log(tau/tau_bar)=(1-bindTax)*(rho_tau*log(tau(-1)/tau_bar)+rho_B*log(B(-1)/B_bar))+bindTax*(log(tau_max));

log(tauSdw) = log(tau_bar) + rho_tau*log(tauSdw(-1)/tau_bar) + rho_B*log(B(-1)/B_bar); 

%Relative:
%log(tau)=(1-bindTax)*(log(tau_bar)+rho_tau*log(tau(-1)/tau_bar)+rho_B*log(B(-1)/Y(-1) * Y_bar/B_bar))+bindTax*(log(tau_max));

%log(tauSdw) = log(tau_bar) + rho_tau*log(tauSdw(-1)/tau_bar) + rho_B*log(B(-1)/Y(-1) * Y_bar/B_bar); 


%[name='8. Government Expenses']
log(G/G_bar)=rho_G*log(G(-1)/G_bar)+rho_Y*log(Y(-1)/Y_bar)+sigma_G*e_G;
log(ExpG/steady_state(G)) = rho_G*log(G/G_bar) ;

%[name='9. Marginal Cost']
%MC=((N^(varphi)*C^(siggma))/(1-tau))*PSI_f*TOT/(A*d_D*Q);
MC=((varphi_N*N^(varphi))/(1-tau))*PSI_f*TOT/(A*K_bar*((G/G_bar)^(gamma_G))*d_D*Q); %GHH

%[name='10. Domestic Production']
N=(Y*d_D)/(A*K_bar*((G/G_bar)^(gamma_G)));

%[name='11. Domestic firms price setting']
x_1_D=(Pi_D/Pi_x_D)*s*(epsilon/(epsilon-1))*x_2_D;

%[name='12. Domestic firms price setting- aux eq 1']
%x_1_D=Gama*(C^(-siggma))*Y+theta_D*betta*((Pi_D(+1)^(epsilon))/Pii(+1))*(Pi_D^((1-epsilon)*delta_D))*x_1_D(+1);
x_1_D=((C+alpha_G*G-(varphi_N*N^(1+varphi))/(1+varphi))^(-siggma))*Y+theta_D*betta*((Pi_D(+1)^(epsilon))/Pii(+1))*(Pi_D^((1-epsilon)*delta_D))*x_1_D(+1); %GHH

%[name='13. Domestic firms price setting- aux eq 2']
%x_2_D=Gama*(C^(-siggma))*MC*Y+theta_D*betta*((Pi_D(+1)^(epsilon))/Pii(+1))*(Pi_D^((-epsilon)*delta_D))*x_2_D(+1);
x_2_D=((C+alpha_G*G-(varphi_N*N^(1+varphi))/(1+varphi))^(-siggma))*MC*Y+theta_D*betta*((Pi_D(+1)^(epsilon))/Pii(+1))*(Pi_D^((-epsilon)*delta_D))*x_2_D(+1);%GHH

%[name='14. Domestic firms Aggregate price']
Pi_D^(1-epsilon)=(1-theta_D)*Pi_x_D^(1-epsilon)+theta_D*Pi_D(-1)^((1-epsilon*delta_D));

%[name='15. Domestic firms price Dispersion']
d_D=(1-theta_D)*(Pi_x_D/Pi_D)^(-epsilon)+theta_D*((Pi_D(-1)^(delta_D)/Pi_D)^(-epsilon))*d_D(-1);

%[name='16. Imports']
M=alppha*(C+G)*(Q/PSI_f)^(-eta);

%[name='17. Importing firms price setting']
x_1_F=(Pi_F/Pi_x_F)*s*(epsilon/(epsilon-1))*x_2_F;

%[name='18. Importing firms price setting- aux eq 1']
%x_1_F=Gama*(C^(-siggma))*M+theta_F*betta*((Pi_F(+1)^(epsilon))/Pii(+1))*(Pi_F^((1-epsilon)*delta_F))*x_1_F(+1);
x_1_F=((C+alpha_G*G(+1)-(varphi_N*N^(1+varphi))/(1+varphi))^(-siggma))*M+theta_F*betta*((Pi_F(+1)^(epsilon))/Pii(+1))*(Pi_F^((1-epsilon)*delta_F))*x_1_F(+1); %GHH

%[name='19. Importing firms price setting- aux eq 2']
%x_2_F=Gama*(C^(-siggma))*M*PSI_f+theta_F*betta*((Pi_F(+1)^(epsilon))/Pii(+1))*(Pi_F^((-epsilon)*delta_F))*x_2_F(+1);
x_2_F=((C+alpha_G*G(+1)-(varphi_N*N^(1+varphi))/(1+varphi))^(-siggma))*M*PSI_f+theta_F*betta*((Pi_F(+1)^(epsilon))/Pii(+1))*(Pi_F^((-epsilon)*delta_F))*x_2_F(+1); %GHH

%[name='20. Importing firms Aggregate price']
Pi_F^(1-epsilon)=(1-theta_F)*Pi_x_F^(1-epsilon)+theta_F*Pi_F(-1)^((1-epsilon*delta_F));

%[name='21. Import firms price Dispersion']
d_F=(1-theta_F)*(Pi_x_F/Pi_F)^(-epsilon)+theta_F*((Pi_F(-1)^(delta_F)/Pi_F)^(-epsilon))*d_F(-1);

%[name='22. Variation- Real exchange rate']
Q=Q(-1)*dS*Pi_x/Pii;

%[name='23. Variation- Terms of trade']
TOT=(TOT(-1))*(Pi_F)/(Pi_D);

%[name='24. Variation- Deviation from LOP']
log(Q)=log(PSI_f)+(1-alppha)*log(TOT);

%[name='25. Inflation relations']
Pii=Pi_D*((TOT/TOT(-1))^alppha);

%[name='26. Market Clearing']
Y=((PSI_f*TOT)^(eta))*(((1-alppha)*(Q^(-eta))*(C+G))+Y_x);

%[name='27. Budget Constraint']
V=((dS/Pii)*R_x(-1)*Phi_v(-1)*V(-1))+((Q*Y)/(PSI_f*TOT))-((Q/(PSI_f*TOT))^(1-eta)*...
(C+G)*(1-alppha+(alppha*PSI_f*TOT^(1-eta))));

%[name='28. Lump-sum transfers']
Z=Z_bar;

%[name='29. Lump-sum taxes']
T_ls=B_bar*((1/R_bar)-(1-delta)/Pi_bar)-G_bar-Z_bar+tau_bar*((C_bar+G_bar)*(1+alppha-alppha*Q_bar));

%Foreign Sector - VAR(1)
%[name='30. External product']
%Y_x=Y_xbar;
log(Y_x/Y_xbar)=a1_yy*log(Y_x(-1)/Y_xbar)+a1_yp*(log(Pi_x(-1)/Pi_xbar))+a1_yi*log(R_x(-1)/R_xbar)+epsilon_yx;

%[name='31. External inflation']
a0_py*log(Y_x/Y_xbar)+log(Pi_x/Pi_xbar)=a1_py*log(Y_x(-1)/Y_xbar)+a1_pp*log(Pi_x(-1)/Pi_xbar)+a1_pi*log(R_x(-1)/R_xbar)+epsilon_pix;
%Pi_x=1;

%[name='32. External interest rate']
a0_iy*log(Y_x/Y_xbar)+a0_ip*log(Pi_x/Pi_xbar)+log(R_x/R_xbar)=a1_iy*log(Y_x(-1)/Y_xbar)+a1_ip*log(Pi_x(-1)/Pi_xbar)+a1_ii*log(R_x(-1)/R_xbar)+epsilon_ix;
%R_x=1/betta;

%[name='33. TFP shock']
log(A/steady_state(A))=rho_A*log(A(-1)/A_bar)+sigma_A*e_A;   
log(ExpA/steady_state(A)) = rho_A*log(A/A_bar) ;


%[name='34. Default']
delta=bindFL*delta_bar;
ExpDelta=delta(+1);

%[name='35. Risk premium shock']
phi=rho_phi*phi(-1)+sigma_phi*e_phi;

%[name='36. Monetary policy shock']
Eps_R=sigma_R*e_R;

%[name='37. External demand shock']
epsilon_yx=sigma_yx*e_yx;

%[name='38. External inflation shock']
epsilon_pix=sigma_pix*e_pix;

%[name='39. International monetary policy shock']
epsilon_ix=sigma_ix*e_ix;


