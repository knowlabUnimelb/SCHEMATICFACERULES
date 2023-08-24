if exist('targ', 'var') ~= 1 || exist('cont', 'var') ~= 1
    fprintf('This code will not work. You need to run computeSIC.m first\n')
    fprintf('(which means that you need to run the analysis code which\n')
    fprintf('calls computeSIC.\n')
    return
end
fig = figure('WindowStyle', 'docked');
subplot(1,2,1)
hold on
e1 = errorbar(1:2, targ(1:2), targerr(1:2), '-k');
set(e1, 'LineWidth', 2)
e2 = errorbar(1:2, targ(3:4), targerr(3:4), '--k');
set(e2, 'LineWidth', 2)

h = plot(1:2, targ(1:2), '-ko', 1:2, targ(3:4), '--ko');
set(gca,'XLim', [.5 2.5], 'XTick', [1 2], 'XTickLabel', {'L', 'H'});
% title(sprintf(['MIC = %4.2f, p =' num2str(p(3))], mic), 'FontSize', 16)
set(h(1), 'MarkerFaceColor', [0 0 0], 'LineWidth', 2, 'MarkerSize',10)
set(h(2), 'MarkerFaceColor', [1 1 1], 'LineWidth', 2, 'MarkerSize',10)
% legend('Low (Top)', 'High (Top)', 'Location', 'NorthEast')
xlabel('Bottom', 'FontSize', 16)
box on

%%
y = anovadata(~ismember(anovadata(:,3), 1:4), [3 5 6 4]);
cont = aggregate(y, 1, 4, [], 1);
contstd = aggregate(y, 1, 4, @std,1); contcnt = aggregate(y, 1, 4, @count,1);
conterr = contstd./sqrt(contcnt);
lowYlim = floor((min([targ; cont] - [targerr; conterr]) - 50)/100) * 100;
highYlim = lowYlim + ceil(max([max([targ; cont] + [targerr; conterr]) - lowYlim + 50, 600])/100) * 100;

set(gca,'YLim', [lowYlim highYlim], 'FontSize', 16)
ylabel('Mean RT (ms)', 'FontSize', 16)

subplot(1,2,2)
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
box on
set(gca,'YLim', [lowYlim highYlim],'FontSize', 16)
xlabel('Interior-Exterior', 'FontSize', 16)
ylabel('Mean RT (ms)', 'FontSize', 16)