close all;
clearvars;
uqlab;
rng(100);

% Input Model
MOpts.mFile = 'model_borehole5D';
myModel = uq_createModel(MOpts);

% Input Marginals
IOpts.Marginals(1).Name = 'rw'; % Radius of the borehole (m)
IOpts.Marginals(1).Type = 'Gaussian';
IOpts.Marginals(1).Parameters = [0.10, 0.0161812];

IOpts.Marginals(2).Name = 'Hu'; % Potentiometric head of the upper aquifer (m)
IOpts.Marginals(2).Type = 'Uniform';
IOpts.Marginals(2).Parameters = [990, 1110];

IOpts.Marginals(3).Name = 'Hl'; % Potentiometric head of the lower aquifer (m)
IOpts.Marginals(3).Type = 'Uniform';
IOpts.Marginals(3).Parameters = [700, 820];

IOpts.Marginals(4).Name = 'L';  % Length of the borehole (m)
IOpts.Marginals(4).Type = 'Uniform';
IOpts.Marginals(4).Parameters = [1120, 1680];

IOpts.Marginals(5).Name = 'Kw'; % Hydraulic conductivity of the borehole (m/yr)
IOpts.Marginals(5).Type = 'Uniform';
IOpts.Marginals(5).Parameters = [9855, 12045];

my_Input = uq_createInput(IOpts);

% Create a PCE Model
PCEOpts.Type = 'Metamodel';
PCEOpts.MetaType = 'PCE';
PCEOpts.FullModel = myModel;
PCEOpts.Degree = 1:10;
PCEOpts.Method = 'LARS';
PCEOpts.ExpDesign.NSamples = 50;
myPCE = uq_createModel(PCEOpts);

% Compute Sobol Indices
%SobolOpts.Type = 'Sensitivity';
%SobolOpts.Method = 'Sobol';
%SobolOpts.Sobol.SampleSize = 1e5;
%PCESobolAnalysis = uq_createAnalysis(SobolOpts);
%PCESobolResults = PCESobolAnalysis.Results;
%STi = PCESobolResults.Total;

% PCE Evaluations for pdf comp

data_pts = load('../valpts_pdf.mat');
X7 = data_pts.S;
X5 = zeros(size(X7,1),5);
X5(:,1) = X7(:,1);
X5(:,2) = X7(:,3);
X5(:,3:5) = X7(:,5:7);
Y_PCE5D = uq_evalModel(X5);
save('y5d.mat','Y_PCE5D');

% Compute UB1 for all variables
X = myPCE.ExpDesign.X;
dim = size(X,2);
nsam = size(X,1);
V = myPCE.PCE.Moments.Var;
dX = zeros(dim,1); pc = zeros(dim,1);
dX(1) = 1e-6; pc(1) = 0.0161812.^2;
for i = 2:dim
    lim = IOpts.Marginals(i).Parameters;
    dX(i) = 1e-5.*(lim(2) - lim(1));
    pc(i) = (1./(pi.^2)).*((lim(2)-lim(1)).^2);
end

%ub1 = dgsm(dim,nsam,X,V,dX,pc);

X7 = load('X7.mat');
Y7 = load('Y7.mat');
X7 = X7.X;
Y7 = Y7.Y;
X5 = zeros(size(X7,1),5);
X5(:,1) = X7(:,1);
X5(:,2) = X7(:,3);
X5(:,3:5) = X7(:,5:7);

errL2 = validate(X5,X7,Y7);

function ub1 = dgsm(dim,nsam,X,V,dX,pc)
G = uq_borehole(X);
Gdx = zeros(nsam,dim);
nu = zeros(dim,1);
ub1 = zeros(dim,1);

for i = 1:dim
    X_new = X;
    X_new(:,i) = X_new(:,i) + dX(i);
    Gdx(:,i) = uq_borehole(X_new);
    dGdx = (Gdx(:,i) - G)./dX(i);
    dGdx2 = dGdx.^2;
    S = sum(dGdx2);
    nu(i) = S./nsam;
    ub1(i) = (pc(i)./V).*nu(i);
end
end

function errv = validate(X5,X7,Y7)
Y_PCE = uq_evalModel(X5);
NR = sum(((Y7-Y_PCE).^2)).^(0.5);
DR = sum((Y7.^2)).^(0.5);
errv = double(NR)./double(DR);
end















