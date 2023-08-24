function [rt, resp] = simmixedSP(n, parms, stoppingrule)

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
parparms = struct('vc', parvc, 'A', Apar, 'bMa', bMa, 's', s, 't0', t0);
[parallelrt, parallelresp] = simparallel(n, parparms, 'st');

serset = [true(nSer,1); false(n-nSer,1)];
parset = [false(nSer,1); true(n-nSer,1)];

rt = [serialrt(serset,1); parallelrt(parset,1)];
resp = [serialresp(serset,1); parallelresp(parset,1)];

rt = rt + t0;            % Add t0