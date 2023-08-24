clear variables
clc
close all
format short g

parallelOn = false; 
whichPool = 'DCS'; % {'DCS', 'local'}

%% Loop through models and subjects
tic
models = {'serialst', 'parallelst', 'serialex', 'parallelex', 'coactive  ', 'mixedSPst', 'mixedSerialC', 'mixedParallelC'};  % Fit models
subjects = [503 504 505 506 601 602 604 606 702 703 704 705 801 802 804 806 807]; % These need to match datafiles in the Data folder
for sidx = 1:numel(subjects)
    for midx = 1:numel(models)
        keep models subjects sidx midx % Clear out variables to avoid any overwriting problems
        subject = subjects(sidx);
       
        %% Model to estimate
        stoppingrule = models{midx}(end-1:end);       % Options: {'st', 'ex', '  '} - Use '  ' for coactive model
        model = models{midx};                         % Options: {'serial', 'parallel', 'coactive', 'mixedSP'}
        fn = sprintf('s%d_%s_t.mat', subject, model); % Set up the save file name
        
        %% Load subject data
        [data, cols] = loadData(subject);
               
        %% Load parameters for inference
        [parms, hyperparms, thetaprior] = loadParmSettings(model);
        names = fieldnames(parms);
        
        %% DE-MCMC sampler settings
        n.parms  = numel(names);
        n.chains = 3 * n.parms; % Rule of thumb: 2 x number of subject level parameters or more

        n.initburn    = 500;    % Number of initial burn in samples
        n.migration   = 200;    % Number of iterations that use migration
        n.postmigburn = 200;    % Number of post-migration burn in samples without migration
        n.mcsamples   = 1000;   % Number of post-burn in samples
        n.burnin    = n.initburn + n.migration + n.postmigburn; % Total burn-in samples
        
        n.migrationStep = 20; % How many trials to run migration
        n.mc        = n.burnin + n.mcsamples; % Total number of samples
        
        %% DE-MCMC sampler parmaeters
        beta    = .001; % Small noise value to add to suggestions
        migprob = .05; % After deterministic migration probabilistically migrate chains with a 1/20 probability
       
        %% Pre-allocate chains
        weight = -Inf * ones(n.chains, n.mc); % Initialize likelihood weight vectors for each chain
        for i = 1:n.parms
            theta.(names{i}) = nan(n.chains, n.mc);
        end
        
        %% Generate some fairly widely-spread random start points
        if exist(fn, 'file') ~= 2
            disp('Initializing Chains')
            if parallelOn
                switch whichPool
                    case 'DCS'
                        %                     p = gcp('nocreate');
                        try parpool('MJSProfile1', n.chains); catch; parpool('local'); end % Open parallel job
                        %                     parpool('MJSProfile1', n.chains);
                    case 'local'
                        try parpool('local'); catch; parpool('local'); end % Open parallel job
                end
            end
            
            for i = 1:n.chains
                while weight(i,1) == -Inf % Ensure valid parms
                    thetaprior = reshape(cell2mat(struct2cell(hyperparms))', 2, n.parms)';
                    for j = 1:n.parms
                        pd = makePriorDistribution(names{j}, thetaprior(j,:)); % Set up the appropriate distribution for each parameter
                        theta.(names{j})(i,1) = random(pd);                    % Sample that parameter
                    end
                    weight(i,1) = logDensLikeLR(getChain(theta, i, 1, 'one'), data, model); % Get likelihood
                end
            end
            starti = 2;
        else
            if parallelOn
                switch whichPool
                    case 'DCS'
                        try parpool('MJSProfile1', n.chains); catch; parpool('local'); end % Open parallel job
                    case 'local'
                        try parpool('local'); catch; parpool('local'); end % Open parallel job
                end
            end
            load(fn, 'i', 'use', 'logtheta', 'theta', 'weight', 't'); % Load existing varaibles
            theta = logtheta;
            starti = i;                                   % Start from last save iteration (may need to update n.mc to take more samples)
        end
            
 
        %% Cycle through remaining samples
        for i = starti:n.mc
            if mod(i, 100) == 0
                fprintf('Iteration %d, Time = %3.2f secs\n', i, t(i-1)); tic; % Display iteration time
            end
            
            use.theta = getChain(theta, i-1, [], 'all');  % Get previous parameters samples
            use.like = weight(:,i-1);                     % Get the previous likelihoods
            
            % Do migration or do cross over
            if i >= n.initburn && i < (n.burnin - n.postmigburn) && mod(i, n.migrationStep) == 0 % Deterministic migration phase
                use = migration(use, data, thetaprior, n, beta, model);         % Migration step
            elseif i > (n.burnin - n.postmigburn) && rand <= migprob                             % Probabilistic migration phase
                use = migration(use, data, thetaprior, n, beta, model);         % Migration step
            else                                                                                 % Cross-over
                inuse = use;
                parfor j = 1:n.chains
                    outuse(j) = crossover(j, inuse, data, thetaprior, n,  beta, model); % Sample new values for each chain
                end
                use = reconstructUse(outuse); clear outuse % Build samples from parallel structure
            end
            weight(:,i) = use.like;                             % Keep new sample target distribution values
            theta = getChain(use.theta, i, [], 'update', theta); % Keep new samples

            t(i) = toc;
        end
        totalTime = sum(t);
        displayTime(totalTime);
        
        delete(gcp('nocreate'))
        
        %% Save sample output
        t(i) = toc;
        logtheta = theta;
        theta = transformSamples(logtheta, data);
        save(fn)
    end
end
