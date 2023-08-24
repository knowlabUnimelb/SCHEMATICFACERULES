if exist('t', 'var') ~= 1 || exist('Rdiff', 'var') ~= 1
    fprintf('This code will not work. You need to run computeSIC.m first\n')
    fprintf('(which means that you need to run the analysis code which\n')
    fprintf('calls computeSIC.\n')
    return
end
fig = figure('WindowStyle', 'docked');
plot(t, zeros(1, length(t)), '-k');
hold on

% Rdiff = cOR_H - cOR_L;
% Rboot = cOR_H_boot - cOR_L_boot;

hc = plot(t, Rdiff); 
set(hc(1), 'Color', 'r', 'LineStyle', '-' , 'LineWidth', 2)
xlabel('t', 'FontSize', 14)
ylabel('R_{DIFF}(t)', 'FontSize', 14)
set(gca,'FontSize', 14, 'XLim', [mint maxt])
% hold on
plot(t, Rdiff+Rboot*sm, 'b', t,  Rdiff-Rboot*sm, 'b');    % plot Bootstrap Confidence Interval (1 std)
set(gca, 'XLim', [300 1200], 'YLim', [-1 1])
title('Resiliency Difference', 'FontSize', 14)