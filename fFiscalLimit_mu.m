function mu = fFiscalLimit_mu(arg1, arg2, varargin)

global fitModel_mu

%diffOrder = 0;

%f_mu  = evalin('base', 'fitModel_mu.PolynomialExpression');
f_mu  = fitModel_mu(arg1,arg2);

%if isa(arg1, 'double')
    mu = f_mu;
end