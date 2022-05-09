# MastersThesis
1) FiscalLimit: Routine to simulate Fiscal Limits and obtain the parameters for endogenous probability

       1.1) FisLim_mcm: Function to Simulate Fiscal Limit following Bi (2014,2016) - Monte Carlo
  
       1.2) parameters_FisLim: Parameters for Fiscal Limit calculation and DSGE model
  
       1.3) FisLim_polyapprox: Approximate Fiscal Limit's mean (mu) and standard deviation by a Polynomial Function
  
       1.4) FisLim_logisticapprox: Fit default probability to a Logistic Function

2) ModelSolver: Endogenous Markov Switching Model Solver - IRFs and simulations
  
       2.1) ModelOpenBi_subsidy.rs: RISE model's script
       2.2) FisLim_steadyStateFile: DSGE steady state file
       
3) Auxiliary codes - Not used yet

       3.1) fFiscalLimit_mu: Function that returns Fiscal Limit's mu for a given state vector
       3.2) fFiscalLimit_std: Function that returns Fiscal Limit's std for a given state vector
       3.3) fFiscalLimit_Prob: Function that returns Fiscal Limit's probability (c.d.f) for a given state vector
