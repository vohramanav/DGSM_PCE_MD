close all
clearvars
uqlab

% Input data
X = dlmread(fullfile('pts_pce19D_20.dat')) ;
Y = dlmread(fullfile('record_id_pce19D_20.dat')) ;
MetaOpts.ExpDesign.X = X;
MetaOpts.ExpDesign.Y = Y;

dim = 19;
L = zeros(1,dim); U = zeros(1,dim); 

% Nominal values for pre-factors a1 to a19 for the 19 reactions

nom = [1.915e14,5.080e04,2.160e08,1.230e04,4.577e19,6.165e15,4.714e18,2.240e22,6.170e19,...
       6.630e13,1.690e14,1.810e13,1.450e16,3.020e12,1.202e17,1.000e13,4.820e13,9.550e06,...
       7.000e12];

L(1,:) = 0.9.*nom(1,:); U(1,:) = 1.1.*nom(1,:);


for ii = 1:dim
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





