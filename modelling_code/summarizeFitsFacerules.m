clear all
clc

subjects = [503 504 505 506 601 602 604 606 702 703 704 705 801 802 804 806 807]; % These need to match datafiles in the Data folder
models = {'serialst', 'parallelst', 'serialex', 'parallelex', 'coactive  ', 'mixedSPst', 'mixedSerialC', 'mixedParallelC'};
subdic1 = nan(numel(subjects),numel(models));
subdic2 = nan(numel(subjects),numel(models));
for sidx = 1:numel(subjects)
    disp(sidx)
    subject = subjects(sidx);
    fitfolder = fullfile(pwd, 'Fits');
    
    dic1 = nan(1,numel(models));
    dic2 = nan(1,numel(models));
    for i = 1:numel(models)
        if exist(fullfile(fitfolder, sprintf('s%d_%s_t.mat', subject, models{i})), 'file') == 2
            load(fullfile(fitfolder, sprintf('s%d_%s_t.mat', subject, models{i})), 'model', 'data', 'theta', 'logtheta', 'weight', 'n')

            n.burnin = n.mc - 750;
            dic1(1,i) = computeDIC(model, data, logtheta, weight(:,n.burnin:end), n.burnin, 1);
            dic2(1,i) = computeDIC(model, data, logtheta, weight(:,n.burnin:end), n.burnin, 2);
            names = fieldnames(theta);
            
        end
    end
    subdic1(sidx,:) = dic1;
    subdic2(sidx,:) = dic2;
end