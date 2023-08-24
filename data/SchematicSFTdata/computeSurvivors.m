function [S, d, acc, t, HI, LO] = computeSurvivors(data, type, t)
% cols of data are item, acc, rt

switch type
    case 'empirical'
        for i = 1:numel(unique(data(:,1)))
            acc{i} = data(data(:,1) == i & data(:,2) == 1, 2);
            d{i}    = data(data(:,1) == i & data(:,2) == 1, 3);
            S(:,i) = 1 - (cumsum(hist(data(data(:,1) == i & data(:,2) == 1, end), t))/nnz(data(data(:,1) == i & data(:,2) == 1, end))');
        end
        HI = nan;
        LO = nan;
    case 'kaplan'
        maxSt = 0;
        for i = 1:numel(unique(data(:,1)))
            d{i} = data(data(:,1) == i, 3);
            [sortd, sidx] = sort(d{i});
            acc{i} = data(data(:,1) == i, 2);
            
            [St{i}, T{i}, Sflo{i}, Sfhi{i}] = ecdf(round(sortd), 'censoring', 1-acc{i}(sidx), 'function' ,'survivor', 'bounds', 'on');
            maxSt = max(maxSt, max(T{i}));
        end
        t = 1:maxSt;
        
        S = nan(maxSt,numel(unique(data(:,1))));
        HI = nan(maxSt,numel(unique(data(:,1))));
        LO = nan(maxSt,numel(unique(data(:,1))));        
        for i = 1:numel(unique(data(:,1)))
            S(1:T{i}(1), i) = 1;
            for j = 2:numel(T{i})
                S(T{i}(j-1):T{i}(j), i) = St{i}(j);
                HI(T{i}(j-1):T{i}(j), i) = Sfhi{i}(j);
                LO(T{i}(j-1):T{i}(j), i) = Sflo{i}(j);                
            end
        end
    case 'nelson'
        % Note that the function "S" below is actually the cumulative
        % Hazard function (i.e., log S)
        maxSt = 0;
        for i = 1:numel(unique(data(:,1)))
            d{i} = data(data(:,1) == i, 3);
            [sortd, sidx] = sort(d{i});
            acc{i} = data(data(:,1) == i, 2);
            
            [St{i}, T{i}, Sflo{i}, Sfhi{i}] = ecdf(round(sortd), 'censoring', 1-acc{i}(sidx), 'function' ,'cumulative hazard', 'bounds', 'on');
            maxSt = max(maxSt, max(T{i}));
        end
        t = 1:maxSt;
        
        S  = nan(maxSt,numel(unique(data(:,1))));
        HI = nan(maxSt,numel(unique(data(:,1))));
        LO = nan(maxSt,numel(unique(data(:,1))));
        for i = 1:numel(unique(data(:,1)))
            S(1:T{i}(1), i) = 1;
            HI(1:T{i}(1), i) = 1;
            LO(1:T{i}(1), i) = 1;
            for j = 2:numel(T{i})
                S(T{i}(j-1):T{i}(j), i) = St{i}(j);
                HI(T{i}(j-1):T{i}(j), i) = Sfhi{i}(j);
                LO(T{i}(j-1):T{i}(j), i) = Sflo{i}(j);
            end
        end        
end