function [ccf, H, t, ccfhi, ccflo] = computeCCF(AH, AL, BH, BL, AHacc, ALacc, BHacc, BLacc, varargin)
% Compute Capacity

% Optional arguments
optargs = {5, 2000, 1000, 'empirical'};
newVals = cellfun(@(x) ~isempty(x), varargin); % skip any new inputs if they are empty
optargs(newVals) = varargin(newVals); % now put these defaults into the valuesToUse cell array, and overwrite the ones specified in varargin.
[mint, maxt, nbootstrap, stype] = optargs{:}; % Place optional args in memorable variable names

%% Compute cdf
t = mint:10:maxt; % #### set t, time vector in msec (MIN : bin size : MAX)

% data = {AH, AL, BH, BL};
% for i = 1:numel(data)
%     cdf(:,i) = cumsum(hist(data{i}, t))/nnz(data{i})';
%     S(:,i) = 1 - cdf(:,i);
% end

data = [[1 * ones(numel(AH),1); 2 * ones(numel(AL), 1); 3 * ones(numel(BH), 1); 4 * ones(numel(BL), 1)],...
    [AHacc; ALacc; BHacc; BLacc],...
    [AH; AL; BH; BL]];
% [S, ~, ~, t] = computeSurvivors(data, 'kaplan', []);
% ccf = (log(S(:,1)) - log(S(:,2))) + (log(S(:,3)) - log(S(:,4)));

[H, ~, ~, t, hi, lo] = computeSurvivors(data, stype, t);
if strcmp(stype, 'empirical')
    S = H;
    H = log(H);
end
H = -H;
hi = -hi;
lo = -lo;
ccf = (H(:,1) - H(:,2)) + (H(:,3) - H(:,4));


switch stype
    case 'kaplan'
        ccfhi = (hi(:,1) - hi(:,2)) + (hi(:,3) - hi(:,4));
        ccflo = (lo(:,1) - lo(:,2)) + (lo(:,3) - lo(:,4));
    case 'empirical'
        [ccf_boot] = f_bootstrap_ccf(S(:,1), S(:,2), S(:,3), S(:,4), t, nbootstrap);    % call the bootstrap function
        ccfhi = ccf + ccf_boot;
        ccflo = ccf - ccf_boot;
end

end


function [std_boot] = f_bootstrap_ccf(AH, AL, BH, BL, t, nbootstrap)

% n iterations of resmapling (with replacement)
% ~~~~~~~~~~~~~~~~~~~~~~~~~~
for n = 1:nbootstrap;
    AH_boot = randsample(AH, length(AH), true);
    AL_boot = randsample(AL, length(AL), true);
    BH_boot = randsample(BH, length(BH), true);
    BL_boot = randsample(BL, length(BL), true);
    
    SaH_boot = 1 - cumsum(hist (AH_boot, t) / nnz(AH_boot));
    SaL_boot = 1 - cumsum(hist (AL_boot, t) / nnz(AL_boot));
    SbH_boot = 1 - cumsum(hist (BH_boot, t) / nnz(BH_boot));
    SbL_boot = 1 - cumsum(hist (BL_boot, t) / nnz(BL_boot));
    
    ccf_boot(:, n) = (log(SaH_boot) - log(SaL_boot)) + (log(SbH_boot) - log(SbL_boot));
end

std_boot = (std(ccf_boot, 0, 2));
end