function [rt, resp] = simlba(n, parms, stoppingrule)
% Generate random samples from the LBA
% This function will only take one item at a time
% Stopping rule is not used for this model

names = fieldnames(parms);
for i = 1:numel(names)
    eval(sprintf('%s = %d;', names{i}, parms.(names{i})));
end
v = [vc, 1 - vc];
b = A + bMa;


%% Simulate model
nchoices=length(v);
% choices = 1:nchoices;

k = A * rand(n, nchoices); % Random starting points

% Sample drift rates
% Since both drift rates are sampled independently, this can result in both
% drifts being negative. 
d = zeros(n,nchoices);
for i = 1:n
    while all(d(i,:) <= 0, 2)
        % Continue to sample random drift rates until at least one is positive
        d(i,:) = v + randn(1, nchoices) * s; 
    end
end

ttf = (b-k)./d;             % Time is distance divided by rate
ttf(ttf < 0) = nan;         % Set any negative times equal to nan
[minrt, resp] = nanmin(ttf, [], 2); % Find the minimum time for an accumulator to finish
rt = minrt + t0;            % Add t0