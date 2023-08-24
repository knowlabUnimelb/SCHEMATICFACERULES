function pdf = getPriorDensity(varname, xval, hyperprior)
% See makePriorDistribution

switch varname
    case {'db1', 'db2', 'sp1', 'sp2','A', 'bMa1', 'bMa2', 's', 't0', 'pX', 'm', 'pSer', 'Aser', 'Apar'}
        pdf = lognormpdf(xval, hyperprior(1), hyperprior(2));
end