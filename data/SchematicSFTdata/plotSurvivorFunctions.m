if exist('t', 'var') ~= 1 || exist('tsf', 'var') ~= 1
    fprintf('This code will not work. You need to run computeSIC.m first\n')
    fprintf('(which means that you need to run the analysis code which\n')
    fprintf('calls computeSIC.\n')
    return
end

fig = figure('WindowStyle', 'docked');
hold on
hs = plot(t, tsf);
set(hs(1), 'Color', 'r', 'LineStyle', '-' , 'LineWidth', 2)
set(hs(2), 'Color', 'r', 'LineStyle', '--', 'LineWidth', 2)
set(hs(3), 'Color', 'b', 'LineStyle', '-' , 'LineWidth', 2)
set(hs(4), 'Color', 'b', 'LineStyle', '--', 'LineWidth', 2)
% legend(hs, 'LL', 'LH', 'HL', 'HH')
xlabel('t', 'FontSize', 16)
ylabel('P (T > t)', 'FontSize', 16)
axis([mint maxt 0 1])
set(gca,'FontSize', 16)
% title('Survivor Functions', 'FontSize', 14)
box on

