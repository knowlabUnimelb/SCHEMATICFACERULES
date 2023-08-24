% v2   - Runs [R] code from Matlab
% v1_3 - Alters set up for running Joe's stats from [R]
% v1_2 - Adds option to use empirical survivor functions for SIC
% v1_0 - Uses censored survivor functions, computes SIC & CCF

% 2017-09-22 DL>> ceSurvivor bug had already been fixed. Did not fix accuracyBug
%   Need to ensure that we report that only min RTs and very long RTs were
%   removed (did not remove outliers by item).
%   Fixed ceOutlier bug but this resulted in no difference due to
%   accuracyBug is:
%        outliers{item} = find(any([itemdata(:,strcmp('rt', cols)) < minrt,...
%               itemdata(:,strcmp('acc', cols)) == 0 & itemdata(:,strcmp('rt', cols)) > means(item,2) + stds(item,2) * 3], 2));
%     ==0 should be == 1
%   Fixed minorMicBug in computation of MIC. Affects mic value shown in figures.

% Set up
clear all
clc
close all force hidden
format short g

stype = 'empirical'; % survivor type: set to 'kaplan' to use censured survivor functions, set to 'empirical' to use empirical histograms
dataPrefix   = '2014_schemecomprules';       % String at the beginning of data file
datalocation = fullfile(pwd, 'rawdata');     % Location folder of datafiles
dataformat   = '%s_s%03d_con0%d_ses0%d.dat'; % Format for datafile name; first string is dataPrefix

modeloutputfolder = fullfile(pwd, 'modeldata'); % Location to save model file
Routputfolder = fullfile(pwd, 'Rdata');         % Location to save data for R analysis

dimensions = {'Eyes', 'Nose-Mouth'}; % Specify descriptive names for your dimensions

si = 1; % which subject to analyse? [Must correspond to the entry in subjectNumber variable]
subjectNumbers   = [503, 504, 505, 506,...
    601, 602, 604, 606,...
    702, 703, 704, 705,...
    801, 802, 804, 806]; % Subject number to analyse (these should match your data file)
conditionNumbers = [5, 5, 5, 5, 6, 6, 6, 6, 7, 7, 7, 7, 8, 8, 8, 8]; % These should match your data file; omit if not needed
sessions = {2:8, 2:8, 2:8, 2:8, 2:8, 2:8, 2:8, 2:8, 2:8, 2:8, 2:8, 2:8, 2:8, 2:8, 2:8, 2:8};  % Specify sessions. You should have one datafile for each session. Can set this to more than one session by stating: 2:5

ploton = true;     % Set to true to display plots
runStats = true;   % Set to true to run statistical analyses
anovaTable = 'on'; % set to 'on' to display the ANOVA table figure

minrt = 200; % Minimum RT cutoff 200ms

% Data file column names (some columns are necessary: 'sub', 'itm', 'rt')
cols = {'sub', 'con', 'rot', 'ses', 'tri', 'itm', 'eyes', 'mouth', 'rsp', 'cat', 'acc', 'rt'};

%% Select from lists defined above
subjectNumber = subjectNumbers(si);
conNumber     = conditionNumbers(si);
session       = sessions{si};

%% Read data
data = [];
for i = 1:numel(session)
    data = [data; dlmread(fullfile(datalocation, sprintf(dataformat, dataPrefix, subjectNumber, conNumber, session(i))))];
end

% Set up block column
nPracTrials = 9;           % Number of practice trials
nTrialsPerSession = 459;   % Number of trials per session
nTrialsPerBlock = [45+9, 45, 45, 45, 45, 45, 45, 45, 45, 45]; % Number of trials in each block

nSessions = numel(session); % Number of sessions
nTrialsPerBlock = repmat(nTrialsPerBlock, 1, nSessions); % Duplicate number of trials per block across sessions
blocks = 1:((size(data, 1)/nSessions - nPracTrials)/45); % Number the blocks
blocks = repmat(blocks, 1, nSessions);                   % Replicate blocks across sessions

% Add in block column to data file
sessionblocks = [];
for i = 1:numel(blocks)
    sessionblocks = [sessionblocks; ones(nTrialsPerBlock(i), 1) * blocks(i)];
end
rdata = [data(:,1:4), sessionblocks,  data(:, 5:end)];
cols = [cols(1:4), 'blk', cols(5:end)]; % Add to col names

%% Throw out practice trials
sessions = unique(rdata(:, strcmp('ses', cols)));
cbd = [];
for i = 1:numel(sessions)
    ddata = rdata(rdata(:,strcmp('ses', cols)) == sessions(i),:);
    if any(ismember(blocks, 1)) % If this block has practice trials, then throw them out
        currentBlockData = ddata(ismember(ddata(:,strcmp('blk', cols)), blocks),:);
        cbd = [cbd; currentBlockData(nPracTrials+1:end,:)];
    else % No practice
        currentBlockData = ddata(ismember(ddata(:,strcmp('blk', cols)), blocks),:);
        cbd = [cbd; currentBlockData];
    end
end
data = cbd;

%% Remove overall long RTs
cutoffPercentiles = repmat([99.9], 1, numel(subjectNumbers)); % Throw out anything greater than this prctile (if 100 then don't throw out anything)
% This is useful for getting rid of extremely long RTs if the P, for instance, took a phonecall or fell asleep

idxTooLong = find(rdata(:,end) > prctile(rdata(:,strcmp(cols, 'rt')), cutoffPercentiles(si)));

rdata(idxTooLong, :) = []; % Delete long RTs
fprintf('Number of long RT trials removed = %d\n', numel(idxTooLong))
condition = rdata(1, strcmp('con', cols));

%% Remove timeouts
idx9 = find(rdata(:,strcmp('acc', cols)) == 9); % Remove timeouts
rdata(idx9,:) = []; % Delete timeouts
fprintf('Number of timeouts removed = %d\n', numel(idx9))

rdata(isnan(rdata(:,strcmp('rt', cols))), :) = []; % Delete nans

%% Convert RTs to msecs if not already in secs
if max(rdata(:,end)) < 1000
    rdata(:,end) = rdata(:,end) * 1000;
end

%% Process data file
% Compute accuracy by blocks
if nSessions == 1
    disp('Accuracy by blocks')
    aggregate(cbd(cbd(:,strcmp('ses', cols)) == session, :), find(strcmp('blk', cols)), find(strcmp('acc', cols)))
    disp('Accuracy this session')
    aggregate(cbd(cbd(:,strcmp('ses', cols)) == sessions,:), find(strcmp('ses', cols)), find(strcmp('acc', cols)))
end

%% Compute averages, remove outliers
x = rdata; % Sort the data by session and by tria

means = aggregate(x, strcmp('itm', cols), strcmp('rt', cols), @nanmean);        % Compute the means for each item
stds  = aggregate(x, strcmp('itm', cols), strcmp('rt', cols), @nanstd);  % Compute the stds for each item

data = rdata;

trialnumber =  (1:size(data,1))';
fintrialnumber =[];
trialdata = [];
cedata = [];
for item = 1:9
    correctAndErrorData = data(data(:,strcmp('itm', cols)) == item,:); % Get data for items.
    cetrial    = trialnumber(data(:,strcmp('itm', cols)) == item,1);
    
    itemdata = data(data(:,strcmp('itm', cols)) == item,:);
    trial    = trialnumber(data(:,strcmp('itm', cols)) == item,1);
    
    % Remove errors
    errors{item} = find(itemdata(:,strcmp('acc', cols)) == 0);
    nerrors(item,1) = numel(errors{item});
    itemdata(errors{item},:) = [];
    trial(errors{item},:) = [];
    
    % Remove outlying RTs < minrt (for errors and correct) or > 3 stds + mean (for correct only)
    outliers{item} = find(any([itemdata(:,strcmp('rt', cols)) < minrt,...
        itemdata(:,strcmp('acc', cols)) == 1 &...
        itemdata(:,strcmp('rt', cols)) > means(item,2) + stds(item,2) * 3], 2));
    noutliers(item,1) = numel(outliers{item});
    itemdata(outliers{item}, :) = [];
    trial(outliers{item}, :) = [];
    
    ceoutliers{item} = find(any([correctAndErrorData(:,strcmp('rt', cols)) < minrt,...
        correctAndErrorData(:,strcmp('acc', cols)) == 1 &...
        correctAndErrorData(:,strcmp('rt', cols)) > means(item,2) + stds(item,2) * 3], 2));
    nceoutliers(item,1) = numel(ceoutliers{item});
    
    correctAndErrorData(ismember(cetrial, cetrial(ceoutliers{item})), :) = []; % Remove outlier rows based on their trial nubmer
    cetrial(ceoutliers{item}, :) = [];
    
    ncorrect(item,1) = size(itemdata,1);
    itemdata(:,strcmp('rt', cols)) = itemdata(:,strcmp('rt', cols));
    mrt(item,1) = mean(itemdata(:,strcmp('rt', cols)));
    
    trialdata = [trialdata; itemdata];
    fintrialnumber = [fintrialnumber; trial];
    cedata = [cedata; correctAndErrorData];
    output(item,:) = [prctile(itemdata(:,strcmp('rt', cols)), [10 30 50 70 90]), ncorrect(item,1), noutliers(item,1), nerrors(item,1), mrt(item,1)];
end
save(fullfile(modeloutputfolder, sprintf('s%d_cedata.mat', subjectNumber)), 'cedata') % Save matfile for model fitting

cleandata = sortrows(trialdata, [find(strcmp('ses', cols)), find(strcmp('blk', cols)), find(strcmp('tri', cols))]);
trialdata(:,strcmp(cols, 'rot') | strcmp(cols, 'blk')) = [];

%% Set up code for [R] analysis
rdata = cedata(:,mstrfind(cols, {'sub', 'con',  'rt', 'acc', 'itm'}));
rdata(:,6:7) = zeros(size(rdata, 1),2);
channelCodes = [2 2; 2 1; 1 2; 1 1; 2 -1; 1 -1; -1 2; -1 1; -1 -1];
for i = 1:size(rdata,1)
    rdata(i, 6:7) = channelCodes(rdata(i,5),:);
end
rdata(:,5) = [];

dlmwrite(fullfile(Routputfolder, sprintf('R_analysis_%s_%d.dat', dataPrefix, subjectNumbers(si))), rdata, 'delimiter', '\t')

%% Collate data for ANOVA
cols = {'sub', 'con', 'ses', 'tri', 'itm', 'top', 'bot', 'rsp', 'cat', 'acc', 'rt'}; %% Copied from top of script

anovadata = trialdata(:, [find(strcmp('sub', cols)), find(strcmp('ses', cols)), find(strcmp('itm', cols)), find(strcmp('rt', cols))]);
anovadata(:, 5) = double(ismember(anovadata(:,3), [1 2 3 4]));
anovadata(ismember(anovadata(:,3), [1 3]), 5) = anovadata(ismember(anovadata(:,3), [1 3]), 5) + 1;
anovadata(:, 6) = double(ismember(anovadata(:,3), [1 2 3 4]));
anovadata(ismember(anovadata(:,3), [1 2]), 6) = anovadata(ismember(anovadata(:,3), [1 2]), 6) + 1;

x = anovadata(ismember(anovadata(:,3), 1:4), [3 5 6 4]);
targ = aggregate(x, [2 3], 4, [],1); mic = targ(1) - targ(2) - targ(3) + targ(4);
targstd = aggregate(x, [2 3], 4, @std,1); targcnt = aggregate(x, [2 3], 4, @count,1);
targerr = targstd./sqrt(targcnt);

%% Run ANOVA
if runStats
    [p, t, stats, terms] = anovan(x(:,4), {x(:,2), x(:,3)}, 'varnames', {'Level 1', 'Level 2'}, 'model', 'full', 'display', anovaTable);
    
    % Run sessions anova
    sessionX = anovadata(ismember(anovadata(:,3), 1:4), [2 3 5 6 4]);
    [anovap, t, stats, terms] = anovan(sessionX(:,5), {sessionX(:,1), sessionX(:,3), sessionX(:,4)},...
        'varnames', {'Session', 'Level 1', 'Level 2'}, 'model', 'full', 'display', anovaTable);
    
    itemComparisons = [5 6; 7 8; 5 9; 6 9; 7 9; 8 9];
    for icIdx = 1:size(itemComparisons, 1)
        item1 = anovadata(ismember(anovadata(:,3), itemComparisons(icIdx,1)), 4);
        item2 = anovadata(ismember(anovadata(:,3), itemComparisons(icIdx,2)), 4);
        [h, ttestp, ci, stats]= ttest2(item1, item2);
        if strcmp(anovaTable, 'on')
            fprintf('Contrast category comparison - item %d vs item %d: t(%d) = %6.2f, p = %3.3f\n',  itemComparisons(icIdx, 1), itemComparisons(icIdx, 2), stats.df, stats.tstat, ttestp);
        end
    end
end
cols = {'sub', 'con', 'ses', 'blk', 'tri', 'itm', 'top', 'bot', 'rsp', 'cat', 'acc', 'rt'};

%% Plot Target Category MICs
if ploton
    fig = figure('WindowStyle', 'docked');
    subplot(3,2,1)
    hold on
    e1 = errorbar(1:2, targ(1:2), targerr(1:2), '-k');
    set(e1, 'LineWidth', 2)
    e2 = errorbar(1:2, targ(3:4), targerr(3:4), '--k');
    set(e2, 'LineWidth', 2)
    
    h = plot(1:2, targ(1:2), '-ko', 1:2, targ(3:4), '--ko');
    set(gca,'XLim', [.5 2.5], 'XTick', [1 2], 'XTickLabel', {'L', 'H'});
    %     title(sprintf(['MIC = %4.2f, p =' num2str(p(3))], mic), 'FontSize', 12)
    set(h(1), 'MarkerFaceColor', [0 0 0], 'LineWidth', 2, 'MarkerSize',10)
    set(h(2), 'MarkerFaceColor', [1 1 1], 'LineWidth', 2, 'MarkerSize',10)
    legend('Low (Nose-Mouth)', 'High (Nose-Mouth)', 'Location', 'NorthEast')
    xlabel('Eyes', 'FontSize', 14)
    box on
    
    %%
    y = anovadata(~ismember(anovadata(:,3), 1:4), [3 5 6 4]);
    cont = aggregate(y, 1, 4, [], 1);
    contstd = aggregate(y, 1, 4, @std,1); contcnt = aggregate(y, 1, 4, @count,1);
    conterr = contstd./sqrt(contcnt);
    lowYlim = floor((min([targ; cont] - [targerr; conterr]) - 50)/100) * 100;
    highYlim = lowYlim + ceil(max([max([targ; cont] + [targerr; conterr]) - lowYlim + 50, 600])/100) * 100;
    
    set(gca,'YLim', [lowYlim highYlim], 'FontSize', 12)
    ylabel('Mean RT (ms)', 'FontSize', 14)
    
    subplot(3,2,2)
    hold on
    e3 = errorbar(1, cont(5), conterr(1)); set(e3, 'LineStyle', 'none', 'Color', [0 0 0]);
    set(e3, 'LineWidth', 2)
    e4 = errorbar(2:3, cont([2 1]), conterr([2 1]), '-k');
    set(e4, 'LineWidth', 2)
    e5 = errorbar(2:3, cont([4 3]), conterr([4 3]), '-k');
    set(e5, 'LineWidth', 2)
    
    h2 = plot(1, cont(5), ' sk', 2:3, cont([4 3]), '-ko', 2:3, cont([2 1]), '-kd');
    
    set(gca,'XLim', [.5 3.5], 'XTick', [1 2 3], 'XTickLabel', {'R', 'I', 'E'});
    set(h2(1), 'MarkerFaceColor', [1 1 1], 'LineWidth', 2, 'MarkerSize',10)
    set(h2(2), 'MarkerFaceColor', [1 1 1], 'LineWidth', 2, 'MarkerSize',10)
    set(h2(3), 'MarkerFaceColor', [0 0 0], 'LineWidth', 2, 'MarkerSize',10)
    
    legend(h2([2 3 1]), dimensions{2}, dimensions{1}, 'Redundant', 'Location', 'NorthWest')
    
    box on
    set(gca,'YLim', [lowYlim highYlim],'FontSize', 12)
    xlabel('Interior-Exterior', 'FontSize', 14)
    ylabel('Mean RT (ms)', 'FontSize', 14)
    title('Contrast Category Mean RTs', 'FontSize', 14)
    cols = {'sub', 'con', 'rot', 'ses', 'blk', 'tri', 'itm', 'sat', 'barpos', 'rsp', 'cat', 'acc', 'rt'}; % copied from line 60
    
    %% Run SIC & CCF Analysis
    data = cleandata;
    data(~ismember(data(:,strcmp('itm', cols)), [1 2 3 4 5 6 7 8 9]), :) = [];
    
    %% Estimate CDF for all items
    mint = min([min(cedata(:,strcmp('rt', cols))), 5]);
    maxt = max([max(cedata(:,strcmp('rt', cols)))]) + 300;
    t = mint:10:maxt; % #### set t, time vector in msec (MIN : bin size : MAX)
    
    items = unique(cedata(:,strcmp('itm', cols)));
    [S, d, acc, tsic] = computeSurvivors(cedata(:,mstrfind(cols, {'itm', 'acc', 'rt'})), stype, t);
    
    %  7  | 3    1
    %     |
    %  8  | 4    2
    %     |--------
    %  9    6    5
    
    
    %% Get MICS
    targ = [mean(d{4}) mean(d{2}) mean(d{3}) mean(d{1})]';
    cont = [mean(d{5}) mean(d{6}) mean(d{7}) mean(d{8}) mean(d{9})]';
    targerr = [std(d{4})/sqrt(length(d{4})) std(d{2})/sqrt(length(d{2})) std(d{3})/sqrt(length(d{3})) std(d{1})/sqrt(length(d{1}))]';
    conterr = [std(d{5})/sqrt(length(d{5})) std(d{6})/sqrt(length(d{6})) std(d{7})/sqrt(length(d{7})) std(d{8})/sqrt(length(d{8})) std(d{9})/sqrt(length(d{9}))]';
    
    %% Target SICs
    % Target SIC item codes: LL = 4, LH = 3, HL = 2, HH = 1; Jun - switched LH to 3 and HL to 2 [Before: LH = 2, HL = 3]
    HH = d{1}; HL = d{2}; LH = d{3}; LL = d{4};
    HHacc = acc{1}; HLacc = acc{2}; LHacc = acc{3}; LLacc = acc{4};
    [sic, tcdf, tsf, tsic, sichi, siclo] = computeSIC(LL, LH, HL, HH, LLacc, LHacc, HLacc, HHacc, mint, maxt, [], stype);
    MIC = mean(LL) - mean(LH) - mean(HL) + mean(HH);
    
    % [Rdiff, Rboot, Rcdfs, RSs] = computeResiliency(d{9}, d{6}, d{8}, d{5}, d{7}, mint, maxt);
    AH = d{5}; AL = d{6}; BH = d{7}; BL = d{8};
    AHacc = acc{5}; ALacc = acc{6}; BHacc = acc{7}; BLacc = acc{8};
    [ccf, ccf_H, tccf, ccfhi, ccflo] = computeCCF(AH, AL, BH, BL, AHacc, ALacc, BHacc, BLacc, mint, maxt);
    % ccfMean = nansum(exp(ccf_H));
    ccfMean = cellfun(@mean, d(5:8));
    Mccf = (ccfMean(1) - ccfMean(2)) + (ccfMean(3) - ccfMean(4));
    
    %% Plot Target Category MICs
    % Plot Survivors
    subplot(3,2,3)
    hs = plot(tsic, tsf);
    set(hs(1), 'Color', 'r', 'LineStyle', '-' , 'LineWidth', 2)
    set(hs(2), 'Color', 'r', 'LineStyle', '--', 'LineWidth', 2)
    set(hs(3), 'Color', 'b', 'LineStyle', '-' , 'LineWidth', 2)
    set(hs(4), 'Color', 'b', 'LineStyle', '--', 'LineWidth', 2)
    legend(hs, 'HH', 'HL', 'LH', 'LL') % Jun - switched HH and LL, and HL and LH around.
    xlabel('t', 'FontSize', 14)
    ylabel('P (T > t)', 'FontSize', 14)
    axis([mint maxt 0 1])
    set(gca,'FontSize', 14)
    title('Survivor Functions', 'FontSize', 14)
    
    % Plot SICS
    subplot(3,2,5)
    hsic = plot(tsic, sic);
    hold on
    set(hsic, 'Color', 'r', 'LineStyle', '-' , 'LineWidth', 2)
    
    hsicCI = plot(tsic, sichi, '--b', tsic, siclo, '--b');
    set(hsicCI, 'LineWidth', 1)
    
    xlabel('t', 'FontSize', 14)
    ylabel('SIC(t)', 'FontSize', 14)
    axis tight
    l = line([mint maxt], [0 0]); set(l, 'Color', 'k')
    set(gca,'FontSize', 14)
    title(sprintf('MIC = %4.2f', MIC), 'FontSize', 14)
    
    %      %% Plot conflict survivors
    %     subplot(3,2,4)
    %     hs = plot(tccf, exp(ccf_H));
    %     set(hs(1), 'Color', 'r', 'LineStyle', '-' , 'LineWidth', 2)
    %     set(hs(2), 'Color', 'r', 'LineStyle', '--', 'LineWidth', 2)
    %     set(hs(3), 'Color', 'b', 'LineStyle', '-' , 'LineWidth', 2)
    %     set(hs(4), 'Color', 'b', 'LineStyle', '--', 'LineWidth', 2)
    %     legend(hs, 'BL', 'BH', 'AL', 'AH')
    %     xlabel('t', 'FontSize', 14)
    %     ylabel('P (T > t)', 'FontSize', 14)
    %     axis([mint maxt 0 1])
    %     set(gca,'FontSize', 14)
    %     title('Survivor Functions', 'FontSize', 14)
    %
    %% Plot conflict contrast function
    sm = 2; % 2 * std boot
    subplot(3,2,6)
    plot(t, zeros(1, length(t)), '-k');
    hold on
    
    hc = plot(tccf, ccf);
    set(hc(1), 'Color', 'r', 'LineStyle', '-' , 'LineWidth', 2)
    hold on
    hcCI = plot(tccf, ccfhi, '--b', tccf,  ccflo, '--b');    % plot 95% Confidence Interval
    set(hcCI, 'LineWidth', 1)
    title(sprintf('Mean_{CCF} = %4.2f', Mccf), 'FontSize', 14)
    xlabel('t', 'FontSize', 14)
    ylabel('CCF(t)', 'FontSize', 14)
    set(gca,'FontSize', 12, 'XLim', [mint maxt])
end
disp('Finished')

%% Generate mean table
meanTable = aggregate(cedata(cedata(:,strcmp(cols, 'acc')) == 1,:), mstrfind(cols, {'itm'}), mstrfind(cols, {'rt'}));
erts = aggregate(cedata(cedata(:,strcmp(cols, 'acc')) == 0,:), mstrfind(cols, {'itm'}), mstrfind(cols, {'rt'}));
totalCount = aggregate(cedata, mstrfind(cols, {'itm'}), mstrfind(cols, {'rt'}), @count);
for i = 1:size(erts,1)
    meanTable(erts(i,1),3) = erts(i,2);
end
meanTable(:,4) = nerrors./aggregate(cedata, mstrfind(cols, {'itm'}), mstrfind(cols, {'rt'}), @count, 1);

%% Try to run [R] stats
rOutputFile = sprintf('R_analysis_%s_%03d.out', dataPrefix, subjectNumber);
os = getenv('OS');
if strcmp(os, 'Windows_NT')
    try
        [status, result] = system(sprintf('R CMD BATCH --no-save --no-restore --slave -"%s" -"%d" flogicalRulesStats.R %s', dataPrefix, subjectNumber, rOutputFile));
    catch
        disp('Failed to run [R] stats. Is [R] on the Windows path?')
    end
else
    disp('R stats not enabled for this OS.')
end
