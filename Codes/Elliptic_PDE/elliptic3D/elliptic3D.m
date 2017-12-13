% close all;
% clearvars;
% uqlab;
% rng(100);
% 
% % Input Model
% MOpts.mFile = 'model_elliptic3D';
% myModel = uq_createModel(MOpts);
% 
% % Input Marginals
% IOpts.Marginals(1).Name = 'c';
% IOpts.Marginals(1).Type = 'Uniform';
% IOpts.Marginals(1).Parameters = [0.9, 1.1];
% 
% IOpts.Marginals(2).Name = 'a';
% IOpts.Marginals(2).Type = 'Uniform';
% IOpts.Marginals(2).Parameters = [0.9, 1.1];
% 
% IOpts.Marginals(3).Name = 'b';
% IOpts.Marginals(3).Type = 'Uniform';
% IOpts.Marginals(3).Parameters = [0.9, 1.1];
% 
% my_Input = uq_createInput(IOpts);
% 
% % Create a PCE Model
% PCEOpts.Type = 'Metamodel';
% PCEOpts.MetaType = 'PCE';
% PCEOpts.FullModel = myModel;
% PCEOpts.Degree = 1:10;
% PCEOpts.Method = 'LARS';
% PCEOpts.ExpDesign.NSamples = 20;
% t = cputime;
% myPCE = uq_createModel(PCEOpts);
% cpu_time = cputime - t;

%==================================================================================
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
%==================================================================================

% Comparison of 3D and 4D PDFs of the observable (Umax)

data_pts = load('valpts_pdf.mat');
data_y4D = load('y4d.mat');

vp = data_pts.S;
Y_PCE4D = data_y4D.Y_PCE4D;

Y_PCE3D = uq_evalModel(vp(:,2:4));

step = (max(Y_PCE4D) - min(Y_PCE4D)).*(1e-4);
pts_PCE4D = min(Y_PCE4D):step:max(Y_PCE4D);
step = (max(Y_PCE3D) - min(Y_PCE3D)).*(1e-4);
pts_PCE3D = min(Y_PCE3D):step:max(Y_PCE3D);
[density_PCE4D,xmesh_PCE4D] = ksdensity(Y_PCE4D,pts_PCE4D);
[density_PCE3D,xmesh_PCE3D] = ksdensity(Y_PCE3D,pts_PCE3D);

figure;
hold on;
plot(xmesh_PCE4D,density_PCE4D,'Linewidth',2,'color','b');
plot(xmesh_PCE3D,density_PCE3D,'Linewidth',2,'color','k');
xlabel('$$\mathrm{U_{max}}$$','interpreter','latex');
ylabel('$$\mathrm{PDF}$$','interpreter','latex');
set(gca,'fontsize',14);
set(gca,'TickLabelInterpreter','latex');
set(gcf,'color',[1,1,1]);
leg = legend({'$\mathrm{4D~PCE~(\mathcal{O}(10^{-5}))}$',...
              '$\mathrm{3D~PCE~(\mathcal{O}(10^{-6}))}$'});
set(leg,'Interpreter','latex');
box on;
print -depsc pdf_comp4Dv3D_elliptic.eps





