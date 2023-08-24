function [sic, cdf, S, t, sichi, siclo] = computeSIC(LL, LH, HL, HH, LLacc, LHacc, HLacc, HHacc, varargin)
% LL, LH, HL and HH are a vectors of data, optional arguments are the minimum and maximum
% response times to use in generating the time variable, t

% Optional arguments
optargs = {5, 2000, 1000, 'empirical'};
newVals = cellfun(@(x) ~isempty(x), varargin); % skip any new inputs if they are empty
optargs(newVals) = varargin(newVals); % now put these defaults into the valuesToUse cell array, and overwrite the ones specified in varargin.
[mint, maxt, nbootstrap, stype] = optargs{:}; % Place optional args in memorable variable names

%% Compute cdf
t = mint:10:maxt; % #### set t, time vector in msec (MIN : bin size : MAX)

data = [[ones(numel(HH),1); 2 * ones(numel(HL), 1); 3 * ones(numel(LH), 1); 4 * ones(numel(LL), 1)],...
        [HHacc; HLacc; LHacc; LLacc],...
        [HH; HL; LH; LL]];
    
[S, data, acc, t, hi, lo] = computeSurvivors(data, stype, t);
cdf = 1 - S;

sic  = S(:,1) - S(:,2) - S(:,3) + S(:,4);
switch stype
    case 'kaplan'
sichi  = hi(:,1) - hi(:,2) - hi(:,3) + hi(:,4);
siclo  = lo(:,1) - lo(:,2) - lo(:,3) + lo(:,4);
    case 'empirical'
  [sic_boot] = f_sicbootstrap(HH, HL, LH, LL, t, nbootstrap);    % call the bootstrap function     
  sichi = sic + sic_boot;
  siclo = sic - sic_boot;
end
        
% [sic_boot] = f_sicbootstrap(HH, HL, LH, LL, t, nbootstrap);    % call the bootstrap function
% sic_boot = 0;
%%
function [std_boot] = f_sicbootstrap(HH, HL, LH, LL, t, nbootstrap)
% f_sicbootstrap returns a bootstrapped sample from the 4 data vectors:
% HH, HL, LH, and LL to compute std_boot, which is the std. confidence interval for SIC(t)
% You can choose other C.I. using percentiles (e.g., 90%)
%
% Ami Eidels, AMPC 2011


% n iterations of resmapling
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    for n = 1:nbootstrap;
        HH_boot = randsample(HH, length(HH), true);
        HL_boot = randsample(HL, length(HL), true);
        LH_boot = randsample(LH, length(LH), true);
        LL_boot = randsample(LL, length(LL), true);
                
        Shh_boot = 1 - cumsum( (hist (HH_boot, t) / nnz(HH_boot)) );
        Shl_boot = 1 - cumsum( (hist (HL_boot, t) / nnz(HL_boot)) );
        Slh_boot = 1 - cumsum( (hist (LH_boot, t) / nnz(LH_boot)) );
        Sll_boot = 1 - cumsum( (hist (LL_boot, t) / nnz(LL_boot)) );

        sic_boot (:, n) = Shh_boot - Shl_boot - Slh_boot + Sll_boot;
        
        % I like the computer to display the value of n, so I know
        % how many iterations I have left -- but this slows down!
%         disp(n)
    end

std_boot = (std(sic_boot, 0, 2));