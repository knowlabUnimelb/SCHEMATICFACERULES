function summarizePosteriors(theta, n, names)

table = cell(numel(names), 4);
for i = 1:n.parms
    temp = theta.(names{i})(:,n.burnin:end);
    temp = temp(:);
    
    table{i,1} = names{i};
    table{i,2} = mode(temp);
    table{i,3} = prctile(temp, [2.5]);
    table{i,4} = prctile(temp, [97.5]);
end
disp(table)