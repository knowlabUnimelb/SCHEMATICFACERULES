function [out, lnL] = logDensLikeLR(x, data, model)
% Organize parameters based on the model
% Use holmes PDA method to compute the likelihood for each data point
% Sum these loglikelihoods to compute the overall log likelihood

%% Set up the model 
if strcmp(model, 'coactive  ') 
    flag = true; 
    stoppingrule = '';
    fmodel = @simcoactive;
elseif strcmp(model, 'mixedSPst')
    flag = false;
    stoppingrule = 'st';
    fmodel = @simmixedSP;
elseif strcmp(model, 'mixedSerialC')
    flag = false;
    stoppingrule = 'st';
    fmodel = @simmixedSerialContaminant;
elseif strcmp(model, 'mixedParallelC')
    flag = false;
    stoppingrule = 'st';
    fmodel = @simmixedParallelContaminant;    
else
    flag = false;
    stoppingrule = model(end-1:end);
    eval(sprintf('fmodel = @sim%s;', model(1:end-2)))
end

%% Use GRT to get drift rates
x.db1 = (logit(x.db1, 'inverse') * (data.stimloc(6,1) - data.stimloc(7,1))) + data.stimloc(7,1);
x.db2 = (logit(x.db2, 'inverse') * (data.stimloc(8,2) - data.stimloc(9,2))) + data.stimloc(9,2);
x.sp1 = exp(x.sp1);
x.sp2 = exp(x.sp2);
x.bMa1 = exp(x.bMa1);
x.bMa2 = exp(x.bMa2);
x.s = exp(x.s);
x.t0 = exp(x.t0);

if ~strcmp(model, 'mixedSPst') && ~strcmp(model, 'mixedSerialC') && ~strcmp(model, 'mixedParallelC')
    x.A = exp(x.A);
    vc = grt_dbt(data.stimloc, [x.db1, x.db2], [x.sp1, x.sp2], flag); % Drift rate for the "A" accumulator
    nitems = size(vc,1);
elseif strcmp(model, 'mixedSPst')
    x.m = exp(x.m);
    x.pSer = logit(x.pSer, 'inverse');
    x.Aser = exp(x.Aser);
    x.Apar = exp(x.Apar);
    servc = grt_dbt(data.stimloc, [x.db1, x.db2], [x.sp1, x.sp2], flag); % Drift rate for the "A" accumulator
    parvc = grt_dbt(data.stimloc, [x.db1, x.db2], x.m * [x.sp1, x.sp2], flag); % Drift rate for the "A" accumulator
    nitems = size(servc,1);
elseif strcmp(model, 'mixedSerialC')
    x.pSer = logit(x.pSer, 'inverse');
    x.Aser = exp(x.Aser);
    servc = grt_dbt(data.stimloc, [x.db1, x.db2], [x.sp1, x.sp2], flag); % Drift rate for the "A" accumulator
    nitems = size(servc,1); 
elseif strcmp(model, 'mixedParallelC')
    x.pPar = logit(x.pSer, 'inverse');
    x.Apar = exp(x.Apar);
    parvc = grt_dbt(data.stimloc, [x.db1, x.db2], [x.sp1, x.sp2], flag); % Drift rate for the "A" accumulator
    nitems = size(parvc,1);     
end

%% Organize the parameters for the model and loop through items
for i = 1:nitems
    if strcmp(model, 'mixedSPst')
        parms(i).('servc') = servc(i,:);
        parms(i).('parvc') = parvc(i,:);
        parms(i).('pSer') = x.('pSer');
        parms(i).('Aser') = x.('Aser');
        parms(i).('Apar') = x.('Apar');
        parms(i).('pX') = logit(x.('pX'), 'inverse');
    elseif strcmp(model, 'mixedSerialC')
        parms(i).('servc') = servc(i,:);
        parms(i).('pSer') = x.('pSer');
        parms(i).('Aser') = x.('Aser');
        parms(i).('pX') = logit(x.('pX'), 'inverse');
        parms(i).('maxRT') = max(data.rt);
        parms(i).('pResp') = mean(data.resp(data.item == i) == 1);
    elseif strcmp(model, 'mixedParallelC')
        parms(i).('parvc') = parvc(i,:);
        parms(i).('pPar') = x.('pPar');
        parms(i).('Apar') = x.('Apar');
        parms(i).('maxRT') = max(data.rt);
        parms(i).('pResp') = mean(data.resp(data.item == i) == 1);           
    elseif strcmp(model, 'serialst')
        parms(i).('pX') = logit(x.('pX'), 'inverse');
        parms(i).('vc') = vc(i,:);
        parms(i).('A') = x.('A');
    else
        parms(i).('vc') = vc(i,:);
        parms(i).('A') = x.('A');
    end
    parms(i).('bMa') = [x.('bMa1'), x.('bMa2')];
    parms(i).('s') = x.('s');
    parms(i).('t0') = x.('t0');
end
%% Loop
fit = nan(nitems,1);
parfor i = 1:nitems
    iparms = parms(i);
    
    % Set up the data for the current item
    idata = [data.resp(data.item == i), data.rt(data.item == i)];
    
    % Get likelihood for each data point
    lik = pdaholmes(fmodel, 5e4, idata, iparms, stoppingrule);
    lnL{i} = lik;
    % Sum over all data points
    fit(i) = sum(log(lik));
end
out = sum(fit);