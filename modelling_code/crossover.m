function use = crossover(idx, use, data, thetaprior, n, beta, model)
% idx   = current chain
% use   = structure containing sampled parameters
% data  = structure containing rt and acc data
% hyper = structure containing hyper prior values
% n     = structure containing number of chains, mc samples, and parms
% beta  = small value for uniform random noise

names = fieldnames(use.theta);
theta = getChain(use.theta, idx, [], 'one');

% Current sample weights
use.weight(idx) = use.like(idx) + logDensPrior(theta, thetaprior); % Weight for the current sample is the target probability of that sample (i.e., likelihood x prior)

% DE-Proposals
gamma = 2.38/sqrt(2 * length(names)); % Tuning parameters
chainIdx = 1:n.chains;
chainIdx(chainIdx == idx) = [];
index = datasample(chainIdx, 2, 'Replace', false); % Select two other chains

for i = 1:numel(names)
    theta.(names{i}) = use.theta.(names{i})(idx,:) +... % Current chain
    gamma * (use.theta.(names{i})(index(1)) - use.theta.(names{i})(index(2))) + ... % gamma-scaled difference between two other chains
    unifrnd(-beta, beta);
end

like = logDensLikeLR(theta, data, model);               % New likelihood
newweight = like + logDensPrior(theta, thetaprior); % New weight from new like and new prior

if isnan(newweight);
    newweight = -Inf; % Replace any nan weights with -Inf
end

if rand < (exp(newweight - use.weight(idx))) % rand < new/old, then update theta and like
    for i = 1:numel(names)
        use.theta.(names{i})(idx,:) = theta.(names{i});
    end
    use.like(idx) = like;
else
    theta = getChain(use.theta, idx, [], 'one');
    use.like(idx) = logDensLikeLR(theta, data, model);
end
