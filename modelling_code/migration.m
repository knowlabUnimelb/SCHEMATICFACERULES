function use = migration(use, data, hyper, n, beta, model)

% Migration set
chainIdx = 1:n.chains;
csetsize = round(n.chains/2);
chainSet = datasample(chainIdx, csetsize, 'Replace', false);
propSet  = circshift(chainSet, [0, 1]);

% Get parameter samples
names = fieldnames(use.theta);
chain(n.chains).theta = [];
for idx = chainSet
    for i = 1:numel(names);
        chain(idx).theta.(names{i}) = use.theta.(names{i})(idx);             % Get parms from current chain
        chain(idx).proptheta.(names{i}) = use.theta.(names{i})(idx) + unifrnd(-beta, beta);             % Get parms from current chain
    end
    
    % Current sample weights
    chain(idx).weight = use.like(idx) + logDensPrior(chain(idx).theta, hyper); % Weight for the current sample is the target probability of that sample (i.e., likelihood x prior)
    
    % Proposal samples weights
    chain(idx).newlike =  logDensLikeLR(chain(idx).proptheta, data, model); % Weight for the current sample is the target probability of that sample (i.e., likelihood x prior)
    chain(idx).newweight =  chain(idx).newlike +...
        logDensPrior(chain(idx).proptheta, hyper); % Weight for the current sample is the target probability of that sample (i.e., likelihood x prior)
end

% Migration
for i = 1:csetsize
    if rand < (exp(chain(propSet(i)).newweight - chain(chainSet(i)).weight))
%         fprintf('Migrating...\n')
        for j = 1:numel(names)
            use.theta.(names{j})(chainSet(i)) = chain(propSet(i)).theta.(names{j});
        end
        use.like(chainSet(i)) = chain(propSet(i)).newlike;
    end
end