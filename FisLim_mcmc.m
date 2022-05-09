%--------------------------------------------------------------------------
% Function to Simulate Fiscal Limit following Bi (2014,2016) - Monte Carlo
% Marina Perrupato Mendonça
% Master Thesis's
%--------------------------------------------------------------------------

function[FisLim]=FisLim_mcmc(A0,G0, Uc0, per, initper, sigma_A, sigma_G,...
    rho_A, rho_G,rho_Y,A_bar, G_bar,betta, epsilon, s, tau_max, varphi_N,...
    varphi,siggma,Z_bar,T_LS,gamma_G,alpha_G,K_bar)

% initialize the shock process
shock_eA = norminv(rand(per,1),0,1)*sigma_A;%norminv(rand(per,1),0,1)*sigma_A; randn(per,1)*(sigma_A);
shock_eG = norminv(rand(per,1),0,1)*sigma_G;%norminv(rand(per,1),0,1)*sigma_G; randn(per,1)*(sigma_G);

% initialize the vector vars
A_vec = zeros(per,1);
A_vec(1) = A0;
G_vec = zeros(per,1);
G_vec(1) = G0;
Uc_vec = zeros(per,1);
Uc_vec(1) = Uc0;

Ymax_vec=zeros(per,1);
Ymax_vec(1)=(((epsilon-1)/(epsilon*s))*((1-tau_max)/varphi_N))^(1/varphi)*(A0*K_bar*(G0/G_bar)^(gamma_G))^(1+(1/varphi));
Cmax_vec=zeros(per,1);
Cmax_vec(1)=Ymax_vec(1)-G0;
Nmax_vec=zeros(per,1);
Nmax_vec(1)=(((epsilon-1)/(epsilon*s))*K_bar*((1-tau_max)*(A0/varphi_N)*(G0/G_bar)^(gamma_G)))^(1/varphi);

T_vec = zeros(per,1);
T_vec(1)=tau_max*Ymax_vec(1);
surp_vec = zeros(per,1);
surp_vec(1)=((betta)^1)*Uc_vec(1)*(T_vec(1) - G0 - Z_bar + T_LS);


% t = 2... infty
% each time period (it)   
for it = 2:per
    
    % productivity and gov spending  
    % A = A_bar^(1-rho_A)*A_vec(it-1)^rho_A*exp(shock_eA(it));
    A =exp((1-rho_A)*log(A_bar) + rho_A*log(A_vec(it-1)) + shock_eA(it));
    A_vec(it) = A;
    
    %G = G_bar^(1-rho_G)*G_vec(it-1)^rho_G*exp(shock_eG(it));
    G =exp((1-rho_G)*log(G_bar) + rho_G*log(G_vec(it-1))+ rho_Y*log(Ymax_vec(it-1)) + shock_eG(it));
    G_vec(it) = G;
    
    % main variables to discounted primary surplus
    Ymax_vec(it)=(((epsilon-1)/(epsilon*s))*((1-tau_max)/varphi_N))^(1/varphi)*(A*K_bar*(G/G_bar)^(gamma_G))^(1+(1/varphi));
    Nmax_vec(it)=(((epsilon-1)/(epsilon*s))*K_bar*((1-tau_max)*A/varphi_N)*(G/G_bar)^(gamma_G))^(1/varphi);
    Cmax_vec(it) = Ymax_vec(it)-G;
    T_vec(it) = tau_max*Ymax_vec(it);
    Uc_vec(it) = (Cmax_vec(it)+(alpha_G*G)-(Nmax_vec(it)^(1+varphi))/(1+varphi))^(-siggma);
    
    % compute discounted primary surplus for period it
  surp_vec(it) = ((betta)^it)*Uc_vec(it)*(T_vec(it) - G - Z_bar + T_LS);
end

FisLim=sum(surp_vec(initper:end)*betta^(-initper)/Uc_vec(initper));

%save A_vec.dat A_vec -ascii -double
%save G_vec.dat G_vec -ascii -double
%save T_vec.dat T_vec -ascii -double
%save surp_vec.dat surp_vec -ascii -double
%save Uc_vec.dat Uc_vec -ascii -double
%save Ymax_vec.dat Ymax_vec -ascii -double
%save Cmax_vec.dat Cmax_vec -ascii -double
%save Nmax_vec.dat Nmax_vec -ascii -double







