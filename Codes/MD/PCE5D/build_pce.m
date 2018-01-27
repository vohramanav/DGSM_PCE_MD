close all
clearvars
uqlab

% Input data
X = dlmread(fullfile('params5D_1_79.dat')) ;
Y = dlmread(fullfile('k_1_79.dat')) ;
MetaOpts.ExpDesign.X = X;
MetaOpts.ExpDesign.Y = Y;

% Input Random Parameters
L = [6.3446 0 1.62 18.90 1.08];
U = [7.7545 0.1 1.98 23.10 1.32];

for ii = 1:5
    IOpts.Marginals(ii).Type = 'Uniform';
    IOpts.Marginals(ii).Parameters = [L(ii), U(ii)];
end

my_Input = uq_createInput(IOpts);

% Set-up PCE
MetaOpts.Type = 'Metamodel';
MetaOpts.MetaType = 'PCE';
MetaOpts.Degree = 1:30;
MetaOpts.Method = 'LARS';

myPCE = uq_createModel(MetaOpts);
uq_print(myPCE);
