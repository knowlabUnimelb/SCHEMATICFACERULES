if exist('t', 'var') ~= 1 || exist('sic', 'var') ~= 1
    fprintf('This code will not work. You need to run computeSIC.m first\n')
    fprintf('(which means that you need to run the analysis code which\n')
    fprintf('calls computeSIC.\n')
    return
end
fig = figure('WindowStyle', 'docked');
hsic = plot(t, sic); 
hold on
set(hsic, 'Color', 'r', 'LineStyle', '-' , 'LineWidth', 2)
hsicCI = plot(t, sic + std_boot', '--b', t, sic - std_boot', '--b');
set(hsicCI, 'LineWidth', 1)

xlabel('t', 'FontSize', 14)
ylabel('SIC(t)', 'FontSize', 16)
axis tight
l = line([mint maxt], [0 0]); set(l, 'Color', 'k')
set(gca,'FontSize', 16)
% title('SIC', 'FontSize', 14)