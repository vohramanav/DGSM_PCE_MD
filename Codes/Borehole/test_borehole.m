close all;
clearvars;
uqlab;
rng(100);

% Input Model
MOpts.mFile = 'uq_borehole';
myModel = uq_createModel(MOpts);

% Input Marginals
IOpts.Marginals(1).Name = 'rw'; % Radius of the borehole (m)
IOpts.Marginals(1).Type = 'Gaussian';
IOpts.Marginals(1).Parameters = [0.10, 0.0161812];

IOpts.Marginals(2).Name = 'Tu'; % Transmissivity of the upper aquifer (m^2/yr)
IOpts.Marginals(2).Type = 'Uniform';
IOpts.Marginals(2).Parameters = [63070, 115600];

IOpts.Marginals(3).Name = 'Hu'; % Potentiometric head of the upper aquifer (m)
IOpts.Marginals(3).Type = 'Uniform';
IOpts.Marginals(3).Parameters = [990, 1110];

IOpts.Marginals(4).Name = 'Tl'; % Transmissivity of the lower aquifer (m^2/yr)
IOpts.Marginals(4).Type = 'Uniform';
IOpts.Marginals(4).Parameters = [63.1, 116];

IOpts.Marginals(5).Name = 'Hl'; % Potentiometric head of the lower aquifer (m)
IOpts.Marginals(5).Type = 'Uniform';
IOpts.Marginals(5).Parameters = [700, 820];

IOpts.Marginals(6).Name = 'L';  % Length of the borehole (m)
IOpts.Marginals(6).Type = 'Uniform';
IOpts.Marginals(6).Parameters = [1120, 1680];

IOpts.Marginals(7).Name = 'Kw'; % Hydraulic conductivity of the borehole (m/yr)
IOpts.Marginals(7).Type = 'Uniform';
IOpts.Marginals(7).Parameters = [9855, 12045];

my_Input = uq_createInput(IOpts);

% Create a PCE Model
PCEOpts.Type = 'Metamodel';
PCEOpts.MetaType = 'PCE';
PCEOpts.FullModel = myModel;
PCEOpts.Degree = 1:10;
PCEOpts.Method = 'LARS';
PCEOpts.ExpDesign.NSamples = 80;
myPCE = uq_createModel(PCEOpts);

% Compute Sobol Indices
%SobolOpts.Type = 'Sensitivity';
%SobolOpts.Method = 'Sobol';
%SobolOpts.Sobol.SampleSize = 1e5;
%PCESobolAnalysis = uq_createAnalysis(SobolOpts);
%PCESobolResults = PCESobolAnalysis.Results;
%STi = PCESobolResults.Total;
%
%% Plot Sobol indices
%uq_figure('Position', [50 50 500 400])
%cm = colormap;
%uq_bar(1:7,PCESobolResults.Total, 0.25,'facecolor', cm(1,:), 'edgecolor', 'none');
%ylabel('$$\mathrm{\mathcal{T}(\theta_i)}$$','interpreter','latex');
%xlabel('$$\mathrm{Parameter}$$','interpreter','latex');
%ylim([0 1]);
%xlim([0 8]);
%xticklabels({'$$\mathrm{r_w}$$','$$\mathrm{T_u}$$','$$\mathrm{H_u}$$','$$\mathrm{T_l}$$',...
%                '$$\mathrm{H_l}$$','$$\mathrm{L}$$','$$\mathrm{K_w}$$'});
%set(gca, 'xtick', 1:length(IOpts.Marginals), 'xticklabel',...
%    xticklabels, 'fontsize',14);
%set(gca,'TickLabelInterpreter','latex');
%grid off;
%set(gcf,'color',[1,1,1]);
%box on;
%
%print -depsc sense_borehole.eps

% Compute UB1 for all variables
%X = myPCE.ExpDesign.X;
%dim = size(X,2);
%nsam = size(X,1);
%%V = myPCE.PCE.Moments.Var;
%dX = zeros(dim,1); pc = zeros(dim,1);
%dX(1) = 1e-6; pc(1) = 0.0161812.^2;
%for i = 2:dim
%    lim = IOpts.Marginals(i).Parameters;
%    dX(i) = 1e-5.*(lim(2) - lim(1));
%    pc(i) = (1./(pi.^2)).*((lim(2)-lim(1)).^2);
%end

%ub1 = dgsm(dim,nsam,X,V,dX,pc);
%data = load('val_pts.mat');
%E = validate(data.val_pts);

function ub1 = dgsm(dim,nsam,X,dX,pc)
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

function errv = validate(vp)
Y_PCE = uq_evalModel(vp);
Y_Model = uq_borehole(vp);
NR = sum(((Y_Model-Y_PCE).^2)).^(0.5);
DR = sum((Y_Model.^2)).^(0.5);
errv = double(NR)./double(DR);
end





















