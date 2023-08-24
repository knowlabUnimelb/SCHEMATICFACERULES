function lik = simlik(varargin)
% Generic function to simulate from a given logical rule model

%% Default parameters
parms.vc = .7; 
parms.bMa = .45;
parms.A = .5;
parms.s = .1;
parms.t0 = .29;

dt = .001;
maxt = 10;

%% Variable input parameters
optargs = {@simlba, 30000, [], parms, dt, maxt, .00001};
newVals = cellfun(@(x) ~isempty(x), varargin); % skip any new inputs if they are empty
optargs(newVals) = varargin(newVals); % now put these defaults into the valuesToUse cell array, and overwrite the ones specified in varargin.
[model, n, data, parms, dt, maxt, lowerLikBound] = optargs{:}; % Place optional args in memorable variable names

%% Simulate from model
[rt, resp] = model(n, parms);

%% Find kernel density for each set of simulations
rts1 = rt(resp == 1);
rts2 = rt(resp == 2);

t = .01:dt:maxt;

h1 = getbandwidth(rts1); % Get bandwidths for each alternative (Turner2014, Eq. 14)
h2 = getbandwidth(rts2);

% Each density has be scaled by the frequency of that alternative; i.e.,
% each density is defective but the sum over discrete response alternatives
% equals 1 (Turner2014, p. 9)
if mean(resp==1) > 0
    try
        dens1 = ksdensity(rts1, t, 'kernel', 'epanechnikov', 'bandwidth', h1) * mean(resp == 1);
    catch
        save temp
        dens1 = ksdensity(rts1, t, 'kernel', 'epanechnikov', 'bandwidth', h1) * mean(resp == 1);
    end
else
    dens1 = zeros(1,numel(t));
end

if mean(resp==2) > 0
    dens2 = ksdensity(rts2, t, 'kernel', 'epanechnikov', 'bandwidth', h2) * mean(resp == 2);
else
    dens2 = zeros(1,numel(t));
end

if isempty(data)            % If no data provided...
    lik = [dens1', dens2']; % Just return the spdf
else                        % Else compute the likelihood of that data
    spdf = [dens1', dens2'];
  
    % Note: failure of this function can be caused by rt data greater than max t
    [~, idx] = ismember(round(data(:,2),3), round(t,3));
    ind = sub2ind(size(spdf), idx, data(:,1));
    lik = spdf(ind);
    lik(lik < lowerLikBound) = lowerLikBound;
end
