% This cell matrix version is hard to use in matlab
clear all
clc

% If it turns out that there is no difference between Set 1 and Set 2 then collapse set 1 and set 2
separateFaceSets = false;

%% Setup stuff
cols = {'Subject', 'Block', 'Number', 'Set', 'Cued', 'Resp', 'Congruent', 'Direction', 'Alignment', 'Study', 'Test', 'Correct', 'Response', 'RT'};
if separateFaceSets
    conditionCols = {'Set', 'Direction', 'Alignment', 'Congruent', 'Resp'};
else
    conditionCols = { 'Direction', 'Alignment', 'Congruent', 'Resp'};
end

%% Read data
datafolder = fullfile(pwd,  'raw_data');
datafiles = dir(fullfile(datafolder, '*.dat'));
data = [];
for i = 1:numel(datafiles)
    subdata = dlmread(fullfile(datafolder, datafiles(i).name));
    subdata(subdata(:, strcmp(cols, 'Response')) == 2,  strcmp(cols, 'Response')) = 0;
    data = [data; subdata];
    
    % Compute d-prime for each subject and each condition
end


hits = aggregate(data(data(:, strcmp(cols, 'Resp')) == 1, :), mstrfind(cols, ['Subject', conditionCols(1:end-1)]), mstrfind(cols, {'Correct'}), @nanmean); % Hit Rate
cr   = aggregate(data(data(:, strcmp(cols, 'Resp')) == 0, :), mstrfind(cols, ['Subject', conditionCols(1:end-1)]), mstrfind(cols, {'Correct'}), @nanmean); % Correct Rejections
counts = aggregate(data(data(:, strcmp(cols, 'Resp')) == 1, :), mstrfind(cols, ['Subject', conditionCols(1:end-1)]), mstrfind(cols, {'Correct'}), @count, 1); % Hit Rate

hitrate = hits(:,end);
hitrate(hitrate == 1) = 1 - 1./sqrt(2 * counts(hitrate == 1));
hitrate(hitrate == 0) = 1./sqrt(2 * counts(hitrate == 0));

farate = 1- cr(:,end);
farate(farate == 1) = 1 - 1./sqrt(2 * counts(farate == 1));
farate(farate == 0) = 1./sqrt(2 * counts(farate == 0));

zHit = norminv(hitrate);
zFA  = norminv(farate);
dprime = zHit - zFA;
c = -.5 * (zHit + zFA);

if separateFaceSets
    means = aggregate([hits(:,1:5), dprime, c], [2 3 4 5], [6 7]);
    stds  = aggregate([hits(:,1:5), dprime, c], [2 3 4 5], [6 7], @std);
else
    means = aggregate([hits(:,[1:4]), dprime, c], [2 3 4], [5 6]);
    stds  = aggregate([hits(:,[1:4]), dprime, c], [2 3 4], [5 6], @std);
end
n = numel(unique(hits(:,1)));

%% Plot set 1 and set 2 separately
plotOrder = [4 3 2 1]; % UA, UM, IA, IM
if separateFaceSets
    for i = 1:2
        setmeans = means(means(:,1) == i, :);
        setstds  = stds(stds(:,1) == i, :);
        cong = setmeans(setmeans(:,4) == 1, :);
        incong = setmeans(setmeans(:,4) == 0, :);
        scong = setstds(setstds(:,4) == 1, :)./sqrt(n);
        sincong = setstds(setstds(:,4) == 0, :)./sqrt(n);
        
        
        subplot(2,2,i)
        h1 = errorbar(1:4, cong(plotOrder, 5), scong(plotOrder, 5), ' ok'); hold on
        h2 = errorbar(1:4, incong(plotOrder, 5), sincong(plotOrder, 5), ' ok');
        x1 = plot(1:4, cong(plotOrder, 5), ' ok', 'MarkerFaceColor', 'k'); hold on
        x2 = plot(1:4, incong(plotOrder, 5), ' ok', 'MarkerFaceColor', 'w');
        set(gca, 'XTick', 1:4, 'XTickLabel', {'UA', 'US', 'IA', 'IS'}, 'FontSize', 10, 'XLim', [.5 4.5])
        xlabel('Condition', 'FontSize', 12)
        ylabel('Sensitivity (d'')', 'FontSize', 12)
        title(sprintf('Set %d', i), 'FontSize', 12)
        legend([x1; x2], 'Congruent', 'Incongruent')
        
        subplot(2,2,i+2)
        h3 = errorbar(1:4, cong(plotOrder, 6), scong(plotOrder, 6), ' ok'); hold on
        h4 = errorbar(1:4, incong(plotOrder, 6), sincong(plotOrder, 6), ' ok');
        x3 = plot(1:4, cong(plotOrder, 6), ' ok', 'MarkerFaceColor', 'k'); hold on
        x4 = plot(1:4, incong(plotOrder, 6), ' ok', 'MarkerFaceColor', 'w');
        line([0 5], [0 0], 'LineStyle', '--', 'Color', 'k')
        set(gca, 'XTick', 1:4, 'XTickLabel', {'UA', 'US', 'IA', 'IS'}, 'FontSize', 10, 'XLim', [.5 4.5])
        xlabel('Condition', 'FontSize', 12)
        ylabel('Criterion (c)', 'FontSize', 12)
        title(sprintf('Set %d', i), 'FontSize', 12)
        legend([x3; x4], 'Congruent', 'Incongruent')
    end
else
        setmeans = means;
        setstds  = stds;
        cong = setmeans(setmeans(:,3) == 1, :);
        incong = setmeans(setmeans(:,3) == 0, :);
        scong = setstds(setstds(:,3) == 1, :)./sqrt(n);
        sincong = setstds(setstds(:,3) == 0, :)./sqrt(n);
        
        
        subplot(1,2,1)
        h1 = errorbar(1:4, cong(plotOrder, 4), scong(plotOrder, 4), ' ok'); hold on
        h2 = errorbar(1:4, incong(plotOrder, 4), sincong(plotOrder, 4), ' ok');
        x1 = plot(1:4, cong(plotOrder, 4), ' ok', 'MarkerFaceColor', 'k'); hold on
        x2 = plot(1:4, incong(plotOrder, 4), ' ok', 'MarkerFaceColor', 'w');
        set(gca, 'XTick', 1:4, 'XTickLabel', {'UA', 'US', 'IA', 'IS'}, 'FontSize', 10, 'XLim', [.5 4.5])
        xlabel('Condition', 'FontSize', 12)
        ylabel('Sensitivity (d'')', 'FontSize', 12)
        title(sprintf('d-prime'), 'FontSize', 12)
        legend([x1; x2], 'Congruent', 'Incongruent')
        
        subplot(1,2,2)
        h3 = errorbar(1:4, cong(plotOrder, 5), scong(plotOrder, 5), ' ok'); hold on
        h4 = errorbar(1:4, incong(plotOrder, 5), sincong(plotOrder, 5), ' ok');
        x3 = plot(1:4, cong(plotOrder, 5), ' ok', 'MarkerFaceColor', 'k'); hold on
        x4 = plot(1:4, incong(plotOrder, 5), ' ok', 'MarkerFaceColor', 'w');
        line([0 5], [0 0], 'LineStyle', '--', 'Color', 'k')
        set(gca, 'XTick', 1:4, 'XTickLabel', {'UA', 'US', 'IA', 'IS'}, 'FontSize', 10, 'XLim', [.5 4.5])
        xlabel('Condition', 'FontSize', 12)
        ylabel('Criterion (c)', 'FontSize', 12)
        title(sprintf('Criterion'), 'FontSize', 12)
        legend([x3; x4], 'Congruent', 'Incongruent')
end

%% ANOVA data
ad0 = [hits(:,1:end-1), dprime, c];

% Three way d-prime, all conditions
[p, table, stats, terms] = anovan(ad0(:,5), mat2cell(ad0(:,1:4), size(ad0,1), [1 1 1 1]),...
    'display', 'on',...
    'model', 'full',...
    'random', 1,...
    'varnames', {'Sub', 'Direction', 'Alignment', 'Congruent'});
effects = {'Direction', 'Alignment', 'Congruent', 'Direction*Alignment',...
    'Direction*Congruent', 'Alignment*Congruent', 'Direction*Alignment*Congruent'};
for i = 1:numel(effects)
    line0{i,1} = getANOVAline(table, effects{i}, 'Sub');
end

% Two way d-prime, UA vs US
ad1 = ad0(ad0(:,2) == 1,:);
[p, table, stats, terms] = anovan(ad1(:,5), mat2cell(ad1(:,[1 3:4]), size(ad1,1), [1 1 1]),...
    'display', 'on',...
    'model', 'full',...
    'random', 1,...
    'varnames', {'Sub', 'Alignment', 'Congruent'});

effects = { 'Alignment', 'Congruent', 'Alignment*Congruent'};
for i = 1:numel(effects)
    line1{i,1} = getANOVAline(table, effects{i}, 'Sub');
end
       

% Two way d-prime, IA vs IS
ad2 = ad0(ad0(:,2) == 0,:);
[p, table, stats, terms] = anovan(ad2(:,5), mat2cell(ad2(:,[1 3:4]), size(ad2,1), [1 1 1]),...
    'display', 'on',...
    'model', 'full',...
    'random', 1,...
    'varnames', {'Sub', 'Alignment', 'Congruent'});
effects = { 'Alignment', 'Congruent', 'Alignment*Congruent'};
for i = 1:numel(effects)
    line2{i,1} = getANOVAline(table, effects{i}, 'Sub');
end

% Two way d-prime, UA vs IA
ad3 = ad0((ad0(:,2) == 1 & ad0(:,3) == 1) | (ad0(:,2) == 0 & ad0(:,3) == 1) ,:);
[p, table, stats, terms] = anovan(ad3(:,5), mat2cell(ad3(:,[1:2 4]), size(ad3,1), [1 1 1]),...
    'display', 'on',...
    'model', 'full',...
    'random', 1,...
    'varnames', {'Sub', 'Direction', 'Congruent'});
effects = { 'Direction', 'Congruent', 'Direction*Congruent'};
for i = 1:numel(effects)
    line3{i,1} = getANOVAline(table, effects{i}, 'Sub');
end

% Two way d-prime, US vs IS
ad4 = ad0((ad0(:,2) == 1 & ad0(:,3) == 0) | (ad0(:,2) == 0 & ad0(:,3) == 0) ,:);
[p, table, stats, terms] = anovan(ad4(:,5), mat2cell(ad4(:,[1:2 4]), size(ad4,1), [1 1 1]),...
    'display', 'on',...
    'model', 'full',...
    'random', 1,...
    'varnames', {'Sub', 'Direction', 'Congruent'});
effects = { 'Direction', 'Congruent', 'Direction*Congruent'};
for i = 1:numel(effects)
    line4{i,1} = getANOVAline(table, effects{i}, 'Sub');
end


% Two way d-prime, UA vs IS
ad5 = ad0((ad0(:,2) == 1 & ad0(:,3) == 1) | (ad0(:,2) == 0 & ad0(:,3) == 0) ,:);
[p, table, stats, terms] = anovan(ad5(:,5), mat2cell(ad5(:,[1:2 4]), size(ad5,1), [1 1 1]),...
    'display', 'on',...
    'model', 'full',...
    'random', 1,...
    'varnames', {'Sub', 'Condition', 'Congruent'});
effects = { 'Condition', 'Congruent', 'Condition*Congruent'};
for i = 1:numel(effects)
    line5{i,1} = getANOVAline(table, effects{i}, 'Sub');
end

% Two way d-prime, US vs IA
ad6 = ad0((ad0(:,2) == 1 & ad0(:,3) == 0) | (ad0(:,2) == 0 & ad0(:,3) == 1) ,:);
[p, table, stats, terms] = anovan(ad6(:,5), mat2cell(ad6(:,[1:2 4]), size(ad6,1), [1 1 1]),...
    'display', 'on',...
    'model', 'full',...
    'random', 1,...
    'varnames', {'Sub', 'Condition', 'Congruent'});
effects = { 'Condition', 'Congruent', 'Condition*Congruent'};
for i = 1:numel(effects)
    line6{i,1} = getANOVAline(table, effects{i}, 'Sub');
end

%% ANOVA data
ad0 = [hits(:,1:end-1), dprime, c];

% Three way d-prime, all conditions
[p, table, stats, terms] = anovan(ad0(:,6), mat2cell(ad0(:,1:4), size(ad0,1), [1 1 1 1]),...
    'display', 'on',...
    'model', 'full',...
    'random', 1,...
    'varnames', {'Sub', 'Direction', 'Alignment', 'Congruent'});
effects = {'Direction', 'Alignment', 'Congruent', 'Direction*Alignment',...
    'Direction*Congruent', 'Alignment*Congruent', 'Direction*Alignment*Congruent'};
for i = 1:numel(effects)
    cline0{i,1} = getANOVAline(table, effects{i}, 'Sub');
end

% Two way d-prime, UA vs US
ad1 = ad0(ad0(:,2) == 1,:);
[p, table, stats, terms] = anovan(ad1(:,6), mat2cell(ad1(:,[1 3:4]), size(ad1,1), [1 1 1]),...
    'display', 'on',...
    'model', 'full',...
    'random', 1,...
    'varnames', {'Sub', 'Alignment', 'Congruent'});

effects = { 'Alignment', 'Congruent', 'Alignment*Congruent'};
for i = 1:numel(effects)
    cline1{i,1} = getANOVAline(table, effects{i}, 'Sub');
end
       

% Two way d-prime, IA vs IS
ad2 = ad0(ad0(:,2) == 0,:);
[p, table, stats, terms] = anovan(ad2(:,6), mat2cell(ad2(:,[1 3:4]), size(ad2,1), [1 1 1]),...
    'display', 'on',...
    'model', 'full',...
    'random', 1,...
    'varnames', {'Sub', 'Alignment', 'Congruent'});
effects = { 'Alignment', 'Congruent', 'Alignment*Congruent'};
for i = 1:numel(effects)
    cline2{i,1} = getANOVAline(table, effects{i}, 'Sub');
end

% Two way d-prime, UA vs IA
ad3 = ad0((ad0(:,2) == 1 & ad0(:,3) == 1) | (ad0(:,2) == 0 & ad0(:,3) == 1) ,:);
[p, table, stats, terms] = anovan(ad3(:,6), mat2cell(ad3(:,[1:2 4]), size(ad3,1), [1 1 1]),...
    'display', 'on',...
    'model', 'full',...
    'random', 1,...
    'varnames', {'Sub', 'Direction', 'Congruent'});
effects = { 'Direction', 'Congruent', 'Direction*Congruent'};
for i = 1:numel(effects)
    cline3{i,1} = getANOVAline(table, effects{i}, 'Sub');
end

% Two way d-prime, US vs IS
ad4 = ad0((ad0(:,2) == 1 & ad0(:,3) == 0) | (ad0(:,2) == 0 & ad0(:,3) == 0) ,:);
[p, table, stats, terms] = anovan(ad4(:,6), mat2cell(ad4(:,[1:2 4]), size(ad4,1), [1 1 1]),...
    'display', 'on',...
    'model', 'full',...
    'random', 1,...
    'varnames', {'Sub', 'Direction', 'Congruent'});
effects = { 'Direction', 'Congruent', 'Direction*Congruent'};
for i = 1:numel(effects)
    cline4{i,1} = getANOVAline(table, effects{i}, 'Sub');
end


% Two way d-prime, UA vs IS
ad5 = ad0((ad0(:,2) == 1 & ad0(:,3) == 1) | (ad0(:,2) == 0 & ad0(:,3) == 0) ,:);
[p, table, stats, terms] = anovan(ad5(:,6), mat2cell(ad5(:,[1:2 4]), size(ad5,1), [1 1 1]),...
    'display', 'on',...
    'model', 'full',...
    'random', 1,...
    'varnames', {'Sub', 'Condition', 'Congruent'});
effects = { 'Condition', 'Congruent', 'Condition*Congruent'};
for i = 1:numel(effects)
    cline5{i,1} = getANOVAline(table, effects{i}, 'Sub');
end

% Two way d-prime, US vs IA
ad6 = ad0((ad0(:,2) == 1 & ad0(:,3) == 0) | (ad0(:,2) == 0 & ad0(:,3) == 1) ,:);
[p, table, stats, terms] = anovan(ad6(:,6), mat2cell(ad6(:,[1:2 4]), size(ad6,1), [1 1 1]),...
    'display', 'on',...
    'model', 'full',...
    'random', 1,...
    'varnames', {'Sub', 'Condition', 'Congruent'});
effects = { 'Condition', 'Congruent', 'Condition*Congruent'};
for i = 1:numel(effects)
    cline6{i,1} = getANOVAline(table, effects{i}, 'Sub');
end