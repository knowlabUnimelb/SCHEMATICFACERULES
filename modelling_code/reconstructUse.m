function use = reconstructUse(outuse)

use.like = diag([outuse.like]);

names = fieldnames(outuse(1).theta);
for j = 1:numel(names)
    temp = [];
    for i = 1:size(outuse, 2)
        temp = [temp, outuse(i).theta.(names{j})];
    end
    use.theta.(names{j}) = diag(temp);
end