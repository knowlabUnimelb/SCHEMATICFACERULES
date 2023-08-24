function line = getANOVAline(table, effect, errorlab)

cols = table(1,:);

dfb = table{strcmp(table(:,strcmp(cols, 'Source')), effect), strcmp(cols, 'd.f.')};
dfw = table{strcmp(table(:,strcmp(cols, 'Source')), sprintf('%s*%s', errorlab, effect)), strcmp(cols, 'd.f.')};
F   = table{strcmp(table(:,strcmp(cols, 'Source')), effect), strcmp(cols, 'F')};
mse = table{strcmp(table(:,strcmp(cols, 'Source')), sprintf('%s*%s', errorlab, effect)), strcmp(cols, 'Mean Sq.')};
p   = table{strcmp(table(:,strcmp(cols, 'Source')), effect), strcmp(cols, 'Prob>F')};
line = sprintf('%d\t %d\t %3.2f\t %3.2f\t %s\t F(%d, %d) = %3.2f, MSE = %3.2f, p %s', dfb, dfw, mse, F, pstring(p), dfb, dfw, F, mse, pstring(p));