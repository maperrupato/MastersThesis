function std = fFiscalLimit_std(x1,x2, varargin)

global fitModel_std

%f_std  = evalin('base', 'fitModel_std.PolynomialExpression');
f_std  = fitModel_std(x1,x2);

%if isa(x1, 'double')
    std = f_std;
end