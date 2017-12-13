close all;
clearvars;
uqlab;
rng(100);

% Input Model
MOpts.mFile = 'model_oscillator';
myModel = uq_createModel(MOpts);

% Input Marginals
IOpts.Marginals(1).Name = 'm'; % Mass
IOpts.Marginals(1).Type = 'Gaussian';
IOpts.Marginals(1).Parameters = [1.0, 0.05];

IOpts.Marginals(2).Name = 'k1'; % Spring Constant
IOpts.Marginals(2).Type = 'Gaussian';
IOpts.Marginals(2).Parameters = [1.0, 0.1];

IOpts.Marginals(3).Name = 'k2'; % Spring Constant
IOpts.Marginals(3).Type = 'Gaussian';
IOpts.Marginals(3).Parameters = [0.1, 0.01];

IOpts.Marginals(4).Name = 'r'; % Displacement
IOpts.Marginals(4).Type = 'Gaussian';
IOpts.Marginals(4).Parameters = [0.5, 0.05];

IOpts.Marginals(5).Name = 'F'; % Force
IOpts.Marginals(5).Type = 'Gaussian';
IOpts.Marginals(5).Parameters = [1.0, 0.2];

IOpts.Marginals(6).Name = 't1'; % time
IOpts.Marginals(6).Type = 'Gaussian';
IOpts.Marginals(6).Parameters = [1.0, 0.2];

my_Input = uq_createInput(IOpts);

% Create a PCE Model
PCEOpts.Type = 'Metamodel';
PCEOpts.MetaType = 'PCE';
PCEOpts.FullModel = myModel;
PCEOpts.Degree = 1:10;
PCEOpts.Method = 'LARS';
PCEOpts.ExpDesign.NSamples = 100;
myPCE = uq_createModel(PCEOpts);

% Compute Sobol Indices
SobolOpts.Type = 'Sensitivity';
SobolOpts.Method = 'Sobol';
SobolOpts.Sobol.SampleSize = 1e5;
PCESobolAnalysis = uq_createAnalysis(SobolOpts);
PCESobolResults = PCESobolAnalysis.Results;
STi = PCESobolResults.Total;

%%% PCE VERIFICATION
%data = load('val_pts.mat');
%E = validate(data.val_pts);


% Sensitivity Plot

%  Comparison plot between the indices: Total Sobol' indices
% Create a nice colormap
 uq_figure('filename','TotalSobolComparison.fig', 'Position', [50 50 500 400]);
 cm = colormap;
 hold on
 uq_bar(1:6,PCESobolResults.Total, 0.25,'facecolor', cm(1,:), 'edgecolor', 'none');
 ylabel('$$\mathrm{Sobol~Total~Effect~Index}$$','interpreter','latex');
 xlabel('$$\mathrm{Parameter}$$','interpreter','latex');
% ylim([0 1])
% xlim([0 7])
 set(gca, 'xtick', 1:length(IOpts.Marginals), 'xticklabel',...
     PCESobolResults.VariableNames, 'fontsize',14);
 set(gcf,'color',[1,1,1]);
 box on;
 print -depsc sens_oscillator.eps
%uq_legend({'MC-based (10^6 simulations)', ...
%    sprintf('PCE-based (%d simulations)', myPCE.ExpDesign.NSamples),...
%    sprintf('LRA-based (%d simulations)', myLRA.ExpDesign.NSamples)},...
%    'location', 'northeast', 'fontsize',14);

function errv = validate(vp)
Y_PCE = uq_evalModel(vp);
Y_Model = uq_borehole(vp);
NR = sum(((Y_Model-Y_PCE).^2)).^(0.5);
DR = sum((Y_Model.^2)).^(0.5);
errv = double(NR)./double(DR);
end

%%% Use 'ub1_conv.m' for the upper bound analysis

% Compute UB1 for all variables
% X = myPCE.ExpDesign.X;
% dim = size(X,2);
% nsam = size(X,1);
% %V = myPCE.PCE.Moments.Var;
% dX = zeros(dim,1); pc = zeros(dim,1);
% for i = 1:dim
%     mu = IOpts.Marginals(i).Parameters(1);
%     sig = IOpts.Marginals(i).Parameters(2);
%     dX(i) = 1e-5.*mu;
%     pc(i) = sig.^2;
% end

%ub1 = dgsm(dim,nsam,X,dX,pc);


% function ub1 = dgsm(dim,nsam,X,dX,pc)
% G = model_oscillator(X);
% Gdx = zeros(nsam,dim);
% nu = zeros(dim,1);
% ub1 = zeros(dim,1);
% 
% for i = 1:dim
%     X_new = X;
%     X_new(:,i) = X_new(:,i) + dX(i);
%     Gdx(:,i) = model_oscillator(X_new);
%     dGdx = (Gdx(:,i) - G)./dX(i);
%     dGdx2 = dGdx.^2;
%     S = sum(dGdx2);
%     nu(i) = S./nsam;
%     ub1(i) = pc(i).*nu(i);
% end
% ub1(1:dim) = ub1(1:dim)./sum(ub1);
% end
























