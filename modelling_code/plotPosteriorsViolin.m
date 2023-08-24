clear all
clc
close all

plotGRT = false;
plotLBA = false;
plotLBA2 = false;
plotLBA3 = true;

% Omit 107 UA3, 306 IA5
% subjects = [ 101 105 107, 109 202 204 206 210 301 302 303 304 306, 401 402 403 404]; % These need to match datafiles in the Data folder
% xlabels = {'UA1', 'UA2', 'UA3', 'UA4', 'UM1', 'UM2', 'UM3', 'UM4', 'IA1', 'IA2', 'IA3', 'IA4', 'IA5', 'IM1', 'IM2', 'IM3', 'IM4'};
subjects = 503;
xlabels = {'U1'};
subjectCutPoints = [4, 8, 13]; % subjectCutPoints

% models = {'mixedSPst'};
models = {'coactive'}
fitfolder = fullfile(pwd, 'Fits');

figure1 = figure('WindowStyle', 'docked');
figure2 = figure('WindowStyle', 'docked');
figure3 = figure('WindowStyle', 'docked');
figure4 = figure('WindowStyle', 'docked');

for sidx = 1:numel(subjects)
    subject = subjects(sidx);
    load(fullfile(fitfolder, sprintf('s%d_%s_t.mat', subject, models{1})), 'model', 'data', 'theta', 'logtheta', 'weight', 'n')
    n.burnin = n.mc - 750;
    
    names = fieldnames(theta);
    for j = 1:numel(names);
        temp = theta.(names{j})(:,n.burnin:end);
        samples.(names{j}) = temp(:);
    end
    
    %% GRT parameters
    if plotGRT
        figure(figure1)
        plot1 = gca;
        currparm = 'db1';
        bandwidth = getbandwidth(samples.(currparm)(:));
        [h.(currparm),L.(currparm),MX.(currparm),MED.(currparm),bw.(currparm)] =...
            violin(samples.(currparm), 'x', sidx-.25,...
            'facecolor', [0 .5 0], 'edgecolor', 'k',...
            'facealpha', .5, 'mc', '', 'medc', '', 'bw', bandwidth);
        hp.(currparm) = plot(sidx-.25, mean(samples.(currparm)), ' ok', 'MarkerSize', 5, 'MarkerFaceColor', 'k');
        hold on
        
        currparm = 'db2';
        bandwidth = getbandwidth(samples.(currparm)(:));
        [h.(currparm),L.(currparm),MX.(currparm),MED.(currparm),bw.(currparm)] =...
            violin(samples.(currparm), 'x', sidx-.25,...
            'facecolor', [0  0 .5], 'edgecolor', 'k',...
            'facealpha', .5, 'mc', '', 'medc', '', 'bw', bandwidth);
        hp.(currparm) = plot(sidx-.25, mean(samples.(currparm)), ' ok', 'MarkerSize', 5, 'MarkerFaceColor', 'k');
        hold on
        
        currparm = 'sp1';
        bandwidth = getbandwidth(samples.(currparm)(:));
        [h.(currparm),L.(currparm),MX.(currparm),MED.(currparm),bw.(currparm)] =...
            violin(samples.(currparm), 'x', sidx-.25,...
            'facecolor', [.5  0 0], 'edgecolor', 'k',...
            'facealpha', .5, 'mc', '', 'medc', '', 'bw', bandwidth);
        hp.(currparm) = plot(sidx-.25, mean(samples.(currparm)), ' ok', 'MarkerSize', 5, 'MarkerFaceColor', 'k');
        hold on
        
        currparm = 'sp2';
        bandwidth = getbandwidth(samples.(currparm)(:));
        [h.(currparm),L.(currparm),MX.(currparm),MED.(currparm),bw.(currparm)] =...
            violin(samples.(currparm), 'x', sidx-.25,...
            'facecolor', [.5  .5 0], 'edgecolor', 'k',...
            'facealpha', .5, 'mc', '', 'medc', '', 'bw', bandwidth);
        hp.(currparm) = plot(sidx-.25, mean(samples.(currparm)), ' ok', 'MarkerSize', 5, 'MarkerFaceColor', 'k');
        hold on
        
        if sidx < numel(subjects)
            if ~ismember(sidx, subjectCutPoints)
                line([sidx + .5 sidx + .5], [0 2.5], 'LineStyle', '--', 'Color', 'k')
            else
                line([sidx + .5 sidx + .5], [0 2.5], 'LineStyle', '-', 'Color', 'k')
            end
        end
    end
    %% LBA parameters
    if plotLBA
        figure(figure2)
        plot2 = gca;
        currparm = 'Aser';
        bandwidth = getbandwidth(samples.(currparm)(:));
        [h.(currparm),L.(currparm),MX.(currparm),MED.(currparm),bw.(currparm)] =...
            violin(samples.(currparm), 'x', sidx-.25,...
            'facecolor', [0 .5 0], 'edgecolor', 'k',...
            'facealpha', .5, 'mc', '', 'medc', '', 'bw', bandwidth);
        hp.(currparm) = plot(sidx-.25, mean(samples.(currparm)), ' ok', 'MarkerSize', 5, 'MarkerFaceColor', 'k');
        hold on
        
        currparm = 'Apar';
        bandwidth = getbandwidth(samples.(currparm)(:));
        [h.(currparm),L.(currparm),MX.(currparm),MED.(currparm),bw.(currparm)] =...
            violin(samples.(currparm), 'x', sidx-.25,...
            'facecolor', [0  0 .5], 'edgecolor', 'k',...
            'facealpha', .5, 'mc', '', 'medc', '', 'bw', bandwidth);
        hp.(currparm) = plot(sidx-.25, mean(samples.(currparm)), ' ok', 'MarkerSize', 5, 'MarkerFaceColor', 'k');
        hold on
        
        currparm = 'bMa1';
        bandwidth = getbandwidth(samples.(currparm)(:));
        [h.(currparm),L.(currparm),MX.(currparm),MED.(currparm),bw.(currparm)] =...
            violin(samples.(currparm), 'x', sidx+.25,...
            'facecolor', [.5  0 0], 'edgecolor', 'k',...
            'facealpha', .5, 'mc', '', 'medc', '', 'bw', bandwidth);
        hp.(currparm) = plot(sidx+.25, mean(samples.(currparm)), ' ok', 'MarkerSize', 5, 'MarkerFaceColor', 'k');
        hold on
        
        currparm = 'bMa2';
        bandwidth = getbandwidth(samples.(currparm)(:));
        [h.(currparm),L.(currparm),MX.(currparm),MED.(currparm),bw.(currparm)] =...
            violin(samples.(currparm), 'x', sidx+.25,...
            'facecolor', [.5  .5 0], 'edgecolor', 'k',...
            'facealpha', .5, 'mc', '', 'medc', '', 'bw', bandwidth);
        hp.(currparm) = plot(sidx+.25, mean(samples.(currparm)), ' ok', 'MarkerSize', 5, 'MarkerFaceColor', 'k');
        hold on
        
        if sidx < numel(subjects)
            if ~ismember(sidx, subjectCutPoints)
                line([sidx + .5 sidx + .5], [0 2.5], 'LineStyle', '--', 'Color', 'k')
            else
                line([sidx + .5 sidx + .5], [0 2.5], 'LineStyle', '-', 'Color', 'k')
            end
        end
    end
    %% Other LBA parameters
    if plotLBA2
        figure(figure3)
        plot3 = gca;
        currparm = 's';
        bandwidth = getbandwidth(samples.(currparm)(:));
        [h.(currparm),L.(currparm),MX.(currparm),MED.(currparm),bw.(currparm)] =...
            violin(samples.(currparm), 'x', sidx-.25,...
            'facecolor', [0 .5 0], 'edgecolor', 'k',...
            'facealpha', .5, 'mc', '', 'medc', '', 'bw', bandwidth);
        hp.(currparm) = plot(sidx-.25, mean(samples.(currparm)), ' ok', 'MarkerSize', 5, 'MarkerFaceColor', 'k');
        hold on
        
        currparm = 't0';
        bandwidth = getbandwidth(samples.(currparm)(:));
        [h.(currparm),L.(currparm),MX.(currparm),MED.(currparm),bw.(currparm)] =...
            violin(samples.(currparm), 'x', sidx+.25,...
            'facecolor', [0  0 .5], 'edgecolor', 'k',...
            'facealpha', .5, 'mc', '', 'medc', '', 'bw', bandwidth);
        hp.(currparm) = plot(sidx+.25, mean(samples.(currparm)), ' ok', 'MarkerSize', 5, 'MarkerFaceColor', 'k');
        hold on
        
        if sidx < numel(subjects)
            if ~ismember(sidx, subjectCutPoints)
                line([sidx + .5 sidx + .5], [0 2.5], 'LineStyle', '--', 'Color', 'k')
            else
                line([sidx + .5 sidx + .5], [0 2.5], 'LineStyle', '-', 'Color', 'k')
            end
        end
    end
    
    %% Other parameters
    if plotLBA3
        figure(figure4)
        plot4 = gca;
        currparm = 'pX';
        bandwidth = getbandwidth(samples.(currparm)(:));
        [h.(currparm),L.(currparm),MX.(currparm),MED.(currparm),bw.(currparm)] =...
            violin(samples.(currparm), 'x', sidx-.25,...
            'facecolor', [0 .5 0], 'edgecolor', 'k',...
            'facealpha', .5, 'mc', '', 'medc', '', 'bw', bandwidth);
        hp.(currparm) = plot(sidx-.25, mean(samples.(currparm)), ' ok', 'MarkerSize', 5, 'MarkerFaceColor', 'k');
        hold on
        
        currparm = 'm';
        pt = prctile(samples.(currparm)(:), 90:100);
        samples.(currparm)(samples.(currparm) > pt(5)) = [];
        bandwidth = getbandwidth(samples.(currparm)(:));
        [h.(currparm),L.(currparm),MX.(currparm),MED.(currparm),bw.(currparm)] =...
            violin(samples.(currparm), 'x', sidx+.25,...
            'facecolor', [0  0 .5], 'edgecolor', 'k',...
            'facealpha', .5, 'mc', '', 'medc', '', 'bw', bandwidth);
        hp.(currparm) = plot(sidx+.25, mean(samples.(currparm)), ' ok', 'MarkerSize', 5, 'MarkerFaceColor', 'k');
        hold on
        
        if sidx < numel(subjects)
            if ~ismember(sidx, subjectCutPoints)
                line([sidx + .5 sidx + .5], [0 5], 'LineStyle', '--', 'Color', 'k')
            else
                line([sidx + .5 sidx + .5], [0 5], 'LineStyle', '-', 'Color', 'k')
            end
        end
    end
end
if plotGRT
    set(plot1, 'XLim', [0 18], 'YLim', [0 2.5], 'XTick', 1:17, 'XTickLabel',...
        xlabels);
    legend([h.db1, h.db2, h.sp1, h.sp2], 'D_{TOP}', 'D_{BOT}', '\sigma_{Top}', '\sigma_{Bot}')
    xlabel('Subject')
    ylabel('Density')
end
if plotLBA
    set(plot2, 'XLim', [0 18], 'YLim', [0 1], 'XTick', 1:17, 'XTickLabel',...
        xlabels);
    legend([h.Aser, h.Apar, h.bMa1, h.bMa2], 'A_{serial}', 'A_{parallel}', 'T_{A}-A', 'T_{B}-A')
    xlabel('Subject')
    ylabel('Density')
end
if plotLBA2
    set(plot3, 'XLim', [0 18], 'YLim', [0 .5], 'XTick', 1:17, 'XTickLabel',...
        xlabels);
    legend([h.s, h.t0], 'Drift Variability', 'Non-decision Time')
    xlabel('Subject')
    ylabel('Density')
end
if plotLBA3
    set(plot4, 'XLim', [0 18], 'YLim', [0 2], 'XTick', 1:17, 'XTickLabel',...
        xlabels);
    legend([h.pX, h.m], 'pX', 'm')
    xlabel('Subject')
    ylabel('Density')
    line([0 18], [1 1], 'LineStyle', '-', 'Color', 'k')
end
