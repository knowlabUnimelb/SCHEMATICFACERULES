% Load separable data
function [data, cols] = loadData(subject)

load(fullfile('Data', 'SchematicRules', sprintf('s%d_cedata.mat', subject))) % Load data
cols = {'sub', 'con', 'rot', 'ses', 'blk', 'tri', 'itm', 'eyes', 'mouth', 'rsp', 'cat', 'acc', 'rt'};

%% MDS information
if subject > 500 && subject < 600
    stimloc = dlmread(fullfile('MDS', 'SchematicRules_1_constrained.crd')); % Stimulus locations
elseif subject > 600 && subject < 700
    stimloc = dlmread(fullfile('MDS', 'SchematicRules_2_constrained.crd')); % Stimulus locations
elseif subject > 700 && subject < 800
    stimloc = dlmread(fullfile('MDS', 'SchematicRules_3_constrained.crd')); % Stimulus locations
else
    stimloc = dlmread(fullfile('MDS', 'SchematicRules_4_constrained.crd')); % Stimulus locations
end
n.items = size(stimloc,1);

subdata = sortrows(cedata, find(strcmp(cols, 'itm')));
data.resp = cedata(:,strcmp(cols, 'rsp'));
data.resp(data.resp == 2) = 0; data.resp = data.resp + 1;
data.rt   = cedata(:,strcmp(cols, 'rt'))/1000;
data.item = cedata(:,strcmp(cols, 'itm'));
data.stimloc = stimloc;
