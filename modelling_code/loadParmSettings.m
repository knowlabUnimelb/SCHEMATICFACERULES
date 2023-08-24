function [parms, hyperparms, thetaprior] = loadParmSettings(model)

%% GRT parms
%% Decision boundary
% db1 = (stimloc(6,1)-stimloc(7,1))/2 + stimloc(7,1); % Decision boundary on x in actual bounded coordinates
% Transform db1 from bounded at [stimloc(7,1), stimloc(6,1)] to t_db1 [-inf, inf]
% To transform back to the db1: db1 = (logit(t_db1, 'inverse') * (stimloc(6,1) - stimloc(7,1))) + stimloc(7,1)
% t_db1 = logit((db1 - stimloc(7,1))/(stimloc(6,1) - stimloc(7,1)))
db1 = 0; % This is t_db1, but just call it db1, and set it halfway 
db1_mu    = 0; db1_sigma = .5;

% db2 = (stimloc(8,2)-stimloc(9,2))/2 + stimloc(9,2); % Decision boundary on y
% Transform db2 from bounded at [stimloc(9,2), stimloc(8,2)] to t_db2 [-inf, inf]
% To transform back to db2: db2 = (logit(x, 'inverse') * (stimloc(8,2) - stimloc(9,2))) + stimloc(9,2)
% t_db2 = logit((db2 - stimloc(9,2))/(stimloc(8,2) - stimloc(9,2)))
db2 = 0;
db2_mu    = 0; db2_sigma = .5;

%% Perceptual variability
% sp is constrained to be greater than 0
sp1 = log(.2);
sp1_mu = -1.5; sp1_sigma = .2;

sp2 = log(.2);
sp2_mu = -1.5; sp2_sigma = .2;

% 'db1'    'db2'    'sp1'    'sp2'    'A'    'bMa1'    'bMa2'    's'    't0'    'pX'
% mus = {1.2, 1.6, .2, .2, .7, .1, .1, .25, .45, 3}; % db1, db2, sp1, sp2, A, b1, b2, s, t0, pX
% sds = { .5,  .5, .5, .5,  .5,  1,  1, .1, .2, 1};

%% LBA parms
%% Startpoint
A    = log(.35);  
A_mu = log(.35); A_sigma = .2;

%% A boundary minus start point
bMa1 = log(.25); % Threshold minus startpoint
bMa1_mu = log(.25); bMa1_sigma = 1;

%% B boundary minus start point
bMa2 = log(.25); % Threshold minus startpoint
bMa2_mu = log(.25); bMa2_sigma = 1;

%% Drift rate variability
s   = log(.25);  % Drift variability
s_mu = log(.25); s_sigma = .5;

%% Non-decision time
t0    = log(.22);  
t0_mu = log(.22); t0_sigma = .2;

%% pX
pX = logit(.5);
pX_mu = logit(.5); pX_sigma = 2;

%% m - Multiplier on perceptual variability
% With these settings, m can grow quite large
% m = log(1);
% m_mu = log(1); m_sigma = 2; 

m = log(1);
m_mu = log(1); m_sigma = .2; 

%% pSer 
pSer = logit(.5);
pSer_mu = logit(.5); pSer_sigma = 2;

%% Aser
Aser    = log(.35);  
Aser_mu = log(.35); Aser_sigma = .2;

%% Apar
Apar    = log(.35);  
Apar_mu = log(.35); Apar_sigma = .2;


%% Parm structure
% Set up structure which hold the parameters for that model (need
% one section for each model with different parameterization)
parms = struct(...
    'db1', db1,...
    'db2', db2,...
    'sp1', sp1,...
    'sp2', sp2,...
    'A', A,...
    'bMa1', bMa1,...
    'bMa2', bMa2,...
    's', s,...
    't0', t0);

hyperparms = struct(...
    'db1_mu', db1_mu,...
    'db1_sigma', db1_sigma,...
    'db2_mu', db2_mu,...
    'db2_sigma', db2_sigma,...
    'sp1_mu', sp1_mu,...
    'sp1_sigma', sp1_sigma,...
    'sp2_mu', sp2_mu,...
    'sp2_sigma', sp2_sigma,...
    'A_mu', A_mu,...
    'A_sigma', A_sigma,...
    'bMa1_mu', bMa1_mu,...
    'bMa1_sigma', bMa1_sigma,...
    'bMa2_mu', bMa2_mu,...
    'bMa2_sigma', bMa2_sigma,...
    's_mu', s_mu,...
    's_sigma', s_sigma,...
    't0_mu', t0_mu,...
    't0_sigma', t0_sigma);

% Each parameter has two hyperparms
thetaprior = [...
    db1_mu, db1_sigma;
    db2_mu, db2_sigma;
    sp1_mu, sp1_sigma;
    sp2_mu, sp2_sigma;
    A_mu, A_sigma;
    bMa1_mu, bMa1_sigma;
    bMa2_mu, bMa2_sigma;
    s_mu, s_sigma;
    t0_mu, t0_sigma];

names = fieldnames(parms);
switch model
    case 'serialst'
        parms = setfield(parms, 'pX', pX);
        hyperparms = setfield(hyperparms, 'pX_mu', pX_mu);
        hyperparms = setfield(hyperparms, 'pX_sigma', pX_sigma);
        thetaprior = [thetaprior; pX_mu, pX_sigma];
    case 'mixedSPst'
        parms = rmfield(parms, 'A');
        parms = setfield(parms, 'pX', pX);
        parms = setfield(parms, 'm', m);
        parms = setfield(parms, 'pSer', pSer);
        parms = setfield(parms, 'Aser', Aser);
        parms = setfield(parms, 'Apar', Apar);
        hyperparms = setfield(hyperparms, 'pX_mu', pX_mu);
        hyperparms = setfield(hyperparms, 'pX_sigma', pX_sigma);        
        hyperparms = rmfield(hyperparms, 'A_mu');
        hyperparms = rmfield(hyperparms, 'A_sigma');
        hyperparms = setfield(hyperparms, 'm_mu', m_mu);
        hyperparms = setfield(hyperparms, 'm_sigma', m_sigma);
        hyperparms = setfield(hyperparms, 'pSer_mu', pSer_mu);
        hyperparms = setfield(hyperparms, 'pSer_sigma', pSer_sigma);
        hyperparms = setfield(hyperparms, 'Aser_mu', Aser_mu);
        hyperparms = setfield(hyperparms, 'Aser_sigma', Aser_sigma);
        hyperparms = setfield(hyperparms, 'Apar_mu', Apar_mu);
        hyperparms = setfield(hyperparms, 'Apar_sigma', Apar_sigma);
        thetaprior(strcmp(names, 'A'), :) = [];
        thetaprior = [thetaprior; pX_mu, pX_sigma; m_mu, m_sigma; pSer_mu, pSer_sigma; Aser_mu, Aser_sigma; Apar_mu, Apar_sigma];
    case 'mixedSerialC'
        parms = rmfield(parms, 'A');
        parms = setfield(parms, 'pX', pX);
        parms = setfield(parms, 'pSer', pSer);
        parms = setfield(parms, 'Aser', Aser);
        hyperparms = setfield(hyperparms, 'pX_mu', pX_mu);
        hyperparms = setfield(hyperparms, 'pX_sigma', pX_sigma);        
        hyperparms = rmfield(hyperparms, 'A_mu');
        hyperparms = rmfield(hyperparms, 'A_sigma');
        hyperparms = setfield(hyperparms, 'pSer_mu', pSer_mu);
        hyperparms = setfield(hyperparms, 'pSer_sigma', pSer_sigma);
        hyperparms = setfield(hyperparms, 'Aser_mu', Aser_mu);
        hyperparms = setfield(hyperparms, 'Aser_sigma', Aser_sigma);
        thetaprior(strcmp(names, 'A'), :) = [];
        thetaprior = [thetaprior; pX_mu, pX_sigma; pSer_mu, pSer_sigma; Aser_mu, Aser_sigma];
    case 'mixedParallelC'
        parms = rmfield(parms, 'A');
        parms = setfield(parms, 'pSer', pSer);
        parms = setfield(parms, 'Apar', Apar);
        hyperparms = rmfield(hyperparms, 'A_mu');
        hyperparms = rmfield(hyperparms, 'A_sigma');
        hyperparms = setfield(hyperparms, 'pSer_mu', pSer_mu);
        hyperparms = setfield(hyperparms, 'pSer_sigma', pSer_sigma);
        hyperparms = setfield(hyperparms, 'Apar_mu', Apar_mu);
        hyperparms = setfield(hyperparms, 'Apar_sigma', Apar_sigma);
        thetaprior(strcmp(names, 'A'), :) = [];
        thetaprior = [thetaprior; pSer_mu, pSer_sigma; Apar_mu, Apar_sigma];         
end