close all;
clearvars;
uqlab;
rng(100);

% Input Model
MOpts.mFile = 'model_elliptic';
myModel = uq_createModel(MOpts);

% Input Marginals
IOpts.Marginals(1).Name = 'k';
IOpts.Marginals(1).Type = 'Uniform';
IOpts.Marginals(1).Parameters = [0.09, 0.11]; 

IOpts.Marginals(2).Name = 'c';
IOpts.Marginals(2).Type = 'Uniform';
IOpts.Marginals(2).Parameters = [0.9, 1.1];

IOpts.Marginals(3).Name = 'a';
IOpts.Marginals(3).Type = 'Uniform';
IOpts.Marginals(3).Parameters = [0.9, 1.1];

IOpts.Marginals(4).Name = 'b';
IOpts.Marginals(4).Type = 'Uniform';
IOpts.Marginals(4).Parameters = [0.9, 1.1];

my_Input = uq_createInput(IOpts);

% Create a PCE Model
PCEOpts.Type = 'Metamodel';
PCEOpts.MetaType = 'PCE';
PCEOpts.FullModel = myModel;
PCEOpts.Degree = 1:10;
PCEOpts.Method = 'LARS';
PCEOpts.ExpDesign.NSamples = 25;
t = cputime;
myPCE = uq_createModel(PCEOpts);
cpu_time = cputime - t;

% ========================================================================
% Compute Sobol Indices
% SobolOpts.Type = 'Sensitivity';
% SobolOpts.Method = 'Sobol';
% SobolOpts.Sobol.SampleSize = 1e5;
% PCESobolAnalysis = uq_createAnalysis(SobolOpts);
% PCESobolResults = PCESobolAnalysis.Results;
% STi = PCESobolResults.Total;
% 
% % Sensitivity Plot
% 
% %  Comparison plot between the indices: Total Sobol' indices
% % Create a nice colormap
%  uq_figure('filename','TotalSobolComparison.fig', 'Position', [50 50 500 400]);
%  cm = colormap;
%  hold on
%  uq_bar(1:4,PCESobolResults.Total, 0.25,'facecolor', cm(1,:), 'edgecolor', 'none');
%  ylabel('$$\mathrm{Sobol~Total~Effect~Index}$$','interpreter','latex');
%  xlabel('$$\mathrm{Parameter}$$','interpreter','latex');
% % ylim([0 1])
% % xlim([0 7])
% 
%  set(gca, 'xtick', 1:length(IOpts.Marginals), 'xticklabel',...
%      PCESobolResults.VariableNames, 'fontsize',14);
%  set(gca,'TickLabelInterpreter','latex');
%  set(gcf,'color',[1,1,1]);
%  box on;
%  print -depsc sens_elliptic.eps

%========================================================================
% Generate a set of pdf evaluations of the observable for comparison with
% the 3D case

%dim = 4; % dimension
%nsams = 1e6; % number of points at which the pdf needs to be evaluated
%
%S = zeros(nsams,dim);
%for j = 1:dim
%    L = IOpts.Marginals(j).Parameters(1);
%    U = IOpts.Marginals(j).Parameters(2);
%    S(:,j) = unifrnd(L,U,nsams,1);
%end
%
%Y_PCE4D = uq_evalModel(S);
%save('valpts_pdf.mat','S');
%save('y4d.mat','Y_PCE4D');







