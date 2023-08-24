function [rt, resp] = simparallel(n, parms, stoppingrule)
% Generate random samples from the serial LBA
% parms should include v (2 dimension x 2 A/B accumulators)

names = fieldnames(parms);
for i = 1:numel(names)
    if numel(parms.(names{i})) > 1
        dstr = repmat('%d,', 1, numel(parms.(names{i})));
        eval(sprintf(['%s(:,1) = [' dstr(1:end-1) ']'';'], names{i}, parms.(names{i})));
    else
        eval(sprintf('%s(:,1) = [%d];', names{i}, parms.(names{i})));
    end
end
v = [vc, 1 - vc];
b = A + bMa;

%% Simulate model
nchoices=size(v,2);

k1 = A * rand(n, nchoices); % Random starting points
k2 = A * rand(n, nchoices); % Random starting points

% Sample drift rates
% Since both drift rates are sampled independently, this can result in both drifts being negative. 
d1 = ones(n,1) * v(1,:) + randn(n, nchoices) * s;
% while any(all(d1 < 0, 2))
%    % Continue to sample random drift rates until at least one is positive
%     d1(all(d1 < 0, 2), :) = ones(sum(all(d1 < 0, 2)),1) * v(1,:) + randn(sum(all(d1 < 0, 2)), nchoices) * s; 
% end

d2 = ones(n,1) * v(2,:) + randn(n, nchoices) * s;
% while any(all(d2 < 0, 2))
%    % Continue to sample random drift rates until at least one is positive
%     d2(all(d2 < 0, 2), :) = ones(sum(all(d2 < 0, 2)),1) * v(2,:) + randn(sum(all(d2 < 0, 2)), nchoices) * s; 
% end
ttfa1            = (b(1)-k1(:,1))./d1(:,1);          % Time is distance divided by rate
ttfa1(ttfa1 < 0)  = Inf;                 % Set any negative times equal to nan
ttfb1            = (b(2)-k1(:,2))./d1(:,2);          % Time is distance divided by rate
ttfb1(ttfb1 < 0)  = Inf;                 % Set any negative times equal to nan
ttf1 = [ttfa1, ttfb1];
[minrt1, resp1] = nanmin(ttf1, [], 2); % Find the minimum time for an accumulator to finish

ttfa2            = (b(1)-k2(:,1))./d2(:,1);          % Time is distance divided by rate
ttfa2(ttfa2 < 0) = Inf;                 % Set any negative times equal to nan
ttfb2            = (b(2)-k2(:,2))./d2(:,2);          % Time is distance divided by rate
ttfb2(ttfb2 < 0) = Inf;                 % Set any negative times equal to nan

ttf2 = [ttfa2, ttfb2];
[minrt2, resp2] = nanmin(ttf2, [], 2); % Find the minimum time for an accumulator to finish

rt   = nan(n,1);
resp = nan(n,1);
if strcmp(stoppingrule, 'ex')
    resp(resp1 == 1 & resp2 == 1) = 1; % Both accumulators end at A
    rt(resp1 == 1 & resp2 == 1) = max([minrt1(resp1 == 1 & resp2 == 1), minrt2(resp1 == 1 & resp2 == 1)], [], 2); % max 
    
    resp(resp1 == 2 & resp2 == 1) = 2; % First accumulator ends at B
    rt(resp1 == 2 & resp2 == 1) = max([minrt1(resp1 == 2 & resp2 == 1), minrt2(resp1 == 2 & resp2 == 1)], [], 2); % max 
    
    resp(resp1 == 1 & resp2 == 2) = 2; % Second accumulator ends at B
    rt(resp1 == 1 & resp2 == 2) = max([minrt1(resp1 == 1 & resp2 == 2), minrt2(resp1 == 1 & resp2 == 2)], [], 2); % max 
    
    resp(resp1 == 2 & resp2 == 2) = 2; % Both accumulators end at B
    rt(resp1 == 2 & resp2 == 2) = max([minrt1(resp1 == 2 & resp2 == 2), minrt2(resp1 == 2 & resp2 == 2)], [], 2); % max 
else % Self-terminating
    resp(resp1 == 1 & resp2 == 1) = 1; % Both accumulators end at A
    rt(resp1 == 1 & resp2 == 1) = max([minrt1(resp1 == 1 & resp2 == 1), minrt2(resp1 == 1 & resp2 == 1)], [], 2); % sum 
    
    resp(resp1 == 2 & resp2 == 1) = 2; % First accumulator ends at B
    rt(resp1 == 2 & resp2 == 1) = minrt1(resp1 == 2 & resp2 == 1);
    
    resp(resp1 == 1 & resp2 == 2) = 2; % Second accumulator ends at B
    rt(resp1 == 1 & resp2 == 2) = minrt2(resp1 == 1 & resp2 == 2);
    
    resp(resp1 == 2 & resp2 == 2) = 2; % Both accumulators end at B
    rt(resp1 == 2 & resp2 == 2) = min([minrt1(resp1 == 2 & resp2 == 2), minrt2(resp1 == 2 & resp2 == 2)], [], 2);
end
rt = rt + t0;            % Add t0