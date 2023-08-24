function transformedSamples = transformSamples(theta, data)

names = fieldnames(theta);
for i = 1:numel(names)
    switch names{i}
        case 'db1'
            transformedSamples.(names{i}) =...
                (logit(theta.(names{i}), 'inverse') *...
                (data.stimloc(6,1) - data.stimloc(7,1))) +...
                data.stimloc(7,1);
        case 'db2'
            transformedSamples.(names{i}) =...
                (logit(theta.(names{i}), 'inverse') *...
                (data.stimloc(8,2) - data.stimloc(9,2))) +...
                data.stimloc(9,2);
        case {'sp1', 'sp2', 'A', 'bMa1', 'bMa2', 's', 't0', 'm', 'Aser', 'Apar'}
            transformedSamples.(names{i}) = exp(theta.(names{i}));
        case {'pX', 'pSer'}
            transformedSamples.(names{i}) = logit(theta.(names{i}), 'inverse');
    end
       
end