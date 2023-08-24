% Plot one figure with all of the parameters in it and then separate
% figures for each parameter
figure('WindowStyle', 'docked')
n.burnin = n.mc - 750;
[nr, nc] = nsubplots(n.parms * 2);
cnt = 1;
for i = 1:n.parms
    subplot(nr, nc, cnt); cnt = cnt + 1;
%     plot(theta.(names{i})', ' .k')
    plot(theta.(names{i})')
    title(names{i})
    ylims = get(gca, 'YLim');
    line([n.burnin n.burnin], ylims, 'LineStyle', '--')
   
   
    subplot(nr, nc, cnt); cnt = cnt + 1;
    temp = theta.(names{i})(:,n.burnin:end);
    [c, e] = hist(temp(:), 20);
    h = bar(e, c./sum(c));
    set(h, 'FaceColor', 'w')
    title(names{i})
end

%%
% figure
% [nr, nc] = nsubplots(n.parms * 2);
% cnt = 1;
% for i = 1:n.parms
%     figure('WindowStyle', 'docked')
%     subplot(1, 2, 1); cnt = cnt + 1;
% %     subplot(nr, nc, cnt); cnt = cnt + 1;
% 
%     plot(1:n.mcsamples, theta.(names{i})(:,n.burnin+1:end)')
%     title(names{i})
%     
%     subplot(1,2,2);
% %     subplot(nr, nc, cnt); cnt = cnt + 1;
%     temp = theta.(names{i})(:,n.burnin:end);
%     [c, e] = hist(temp(:), 20);
%     h = bar(e, c./sum(c));
%     set(h, 'FaceColor', 'w')
%     
% end

% %% Posterior predictive check
% figure
% [nr, nc] = nsubplots(n.items);
% cnt = 1;
% for k = 1:n.items
%     subplot(nr,nc,k); 
%     n.simsubs = 20;
%     simdata.rt = []; simdata.acc = [];
%     for i = round(linspace(n.burnin, n.mc, n.simsubs))
%         x.('vc') = theta.(sprintf('vc%d', k))(datasample(1:n.chains, 1),i);
%         for j = (n.items+1):n.parms
%             x.(names{j}) = theta.(names{j})(datasample(1:n.chains, 1),i);
%         end
%         [tmp.rt, tmp.resp] = simlba(nd, x);
%         simdata.rt = [simdata.rt; tmp.rt];
%         simdata.acc = [simdata.acc; tmp.resp == 1];
%     end
%     tmp.datart = data.rt;
%     tmp.datart(~data.acc)=-1*tmp.datart(~data.acc);
%     tmp.datart(abs(tmp.datart)>3)=nan; % Just remove some outliers for plotting purposes.
%     tmp.dataacc = data.acc;
%     
%     [c, e] = hist(tmp.datart, 50);
%     bar(e, c./trapz(e,c), 'hist'); % * mean(tmp.dataacc == 1), 'hist');
%     hold on
%     
%     tmp.simrt = simdata.rt;
%     tmp.simrt(~simdata.acc) = -1*tmp.simrt(~simdata.acc);
%     tmp.simrt(abs(tmp.simrt)>3)=nan; % Just remove some outliers for plotting purposes.
%     
%     hbw = getbandwidth(tmp.simrt(~isnan(tmp.simrt)));
%     [dens, xi] = ksdensity(tmp.simrt(~isnan(tmp.simrt)), 'kernel', 'epanechnikov', 'bandwidth', .01); % *...
%     plot(xi, dens, '-r', 'LineWidth', 2)
% end