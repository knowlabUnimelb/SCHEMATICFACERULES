function y = logTruncatedNormalDistribution(x, mu, sigma, a, b)

if (x > a) & (x < b)
    x1 = (x - mu)./sigma;
    f1 = -log(sigma) + logStandardNormalPdf(x1);
    
    if a == -inf
        f2 = 0;
    else
        a1 = (a - mu)./sigma;
        f2 = standardNormalCdf(a1);
    end
    
    if b == inf
        f3 = 1;
    else
        b1 = (b - mu)./sigma;
        f3 = standardNormalCdf(b1);
    end
    
    y = f1 - log(f2 + f3);
else
    y = -Inf;
end

function y = logStandardNormalPdf(x)
y = (-0.5 * x.^2) - .5 * log(2*pi);

function y = standardNormalCdf(x)
% Use the complementary error function, rather than .5*(1+erf(z/sqrt(2))),
% to produce accurate near-zero results for large negative x.
y = (0.5 * erfc(-x ./ sqrt(2)));
