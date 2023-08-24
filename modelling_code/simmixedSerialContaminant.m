function [rt, resp] = simmixedSerialContaminant(n, parms, stoppingrule)

names = fieldnames(parms);
for i = 1:numel(names)
    if numel(parms.(names{i})) > 1
        dstr = repmat('%d,', 1, numel(parms.(names{i})));
        eval(sprintf(['%s(:,1) = [' dstr(1:end-1) ']'';'], names{i}, parms.(names{i})));
    else
        eval(sprintf('%s(:,1) = [%d];', names{i}, parms.(names{i})));
    end
end

nSer = round(n * pSer);

%% Serial Samples
serparms = struct('vc', servc, 'A', Aser, 'bMa', bMa, 's', s, 't0', t0, 'pX', pX);
[serialrt, serialresp] = simserial(n, serparms, 'st');

%% Parallel Samples
% parparms = struct('vc', parvc, 'A', Apar, 'bMa', bMa, 's', s, 't0', t0);
% [contaminantrt, contaminantresp] = simparallel(n, parparms, 'st');

%% Contaminant Samples
contaminantrt = unifrnd(0, maxRT, n, 1);
contaminantresp = (rand(n, 1) < pResp);
contaminantresp(contaminantresp==0) = 2;

serset = [true(nSer,1); false(n-nSer,1)];
cset   = [false(nSer,1); true(n-nSer,1)];

rt = [serialrt(serset,1); contaminantrt(cset,1)];
resp = [serialresp(serset,1); contaminantresp(cset,1)];

rt = rt + t0;            % Add t0