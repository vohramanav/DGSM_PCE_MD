close all
clearvars
uqlab

% Input data
X = dlmread(fullfile('params5D_1_153.dat')) ;
Y = dlmread(fullfile('k_1_153.dat')) ;
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

% Plots
err_samples();

nsams = 35;
qoi = load('k_data35.txt');
qoic = zeros(nsams-1,1); % ignoring the faulty point
qoic(1:23,1) = qoi(1:23,2);
qoic(24:nsams-1,1) = qoi(25:nsams,2);

%E1 = verify_L2(nsams,qoic);
%verify_pdf(qoic);
%gsa();

% Sobol Sensitivity Indices
function sen = gsa()
SobolOpts.Type = 'Sensitivity';
SobolOpts.Method = 'Sobol';
SobolOpts.Sobol.SampleSize = 1e6;
PCESobolAnalysis = uq_createAnalysis(SobolOpts);
PCESobolResults = PCESobolAnalysis.Results;
TO = PCESobolResults.Total;
FO = PCESobolResults.FirstOrder;
SI = zeros(5,2);
SI(:,1) = FO;
SI(:,2) = TO;

figure;
bar(SI);
ylabel('$$\mathrm{Sensitivity}$$','interpreter','latex');
xtickl = ({'$$\mathrm{A}$$','$$\mathrm{q}$$','$$\mathrm{\alpha}$$',...
    '$$\mathrm{\lambda}$$','$$\mathrm{\gamma}$$','interpreter','latex'});
set(gca, 'xtick', 1:5, 'xticklabel',xtickl, 'fontsize',14);
set(gca,'TickLabelInterpreter','latex');
set(gcf,'color',[1,1,1]);
print -depsc PCE5D_gsa.eps
end

% Plot e_loo vs samples
function err = err_samples()
data = load('err_samples.txt');
semilogy(data(:,1),data(:,2),'--o');
xlabel('$$\mathrm{Number~of~Realizations}$$','interpreter','latex','fontsize',20);
ylabel('$$\mathrm{\log_{10}(\epsilon_{LOO})}$$','interpreter','latex','fontsize',20);
xlim([30,160]);
%ylim([1e-3,1e0]);
set(gca,'TickLabelInterpreter','latex','fontsize',18);
set(gcf,'color',[1,1,1]);
grid on;
print -depsc PCE5D_eloo.eps
end

% PCE Verification
function errl = verify_L2(nsams,qoic)
pts = load('params_35.txt');
ptsc = zeros(nsams-1,5); % ignoring faulty point and using 5D germ
ptsc(1:23,1) = pts(1:23,1);
ptsc(1:23,2:5) = pts(1:23,4:7);
ptsc(24:nsams-1,1) = pts(25:nsams,1);
ptsc(24:nsams-1,2:5) = pts(25:nsams,4:7);

Y_PCE = uq_evalModel(ptsc);
Y_Model = qoic;
NR = sum(((Y_Model-Y_PCE).^2)).^(0.5);
DR = sum((Y_Model.^2)).^(0.5);
errl = double(NR)./double(DR);
end

function errp = verify_pdf(qoic)
pts = load('samples1e6.txt');
ptsc = zeros(1e6,5);
ptsc(:,1) = pts(:,1);
ptsc(:,2:5) = pts(:,4:7);
Y_PCE = uq_evalModel(ptsc);
step = (max(Y_PCE) - min(Y_PCE)).*(1e-4);
pts_PCE = min(Y_PCE):step:max(Y_PCE);
[density_PCE,xmesh_PCE] = ksdensity(Y_PCE,pts_PCE);

figure;
hold on;
nbins = 15;
histogram(qoic,nbins,'Normalization','pdf','FaceAlpha',0.2);
plot(xmesh_PCE,density_PCE,'Linewidth',2,'color','b');
%plot(xmesh_Model,density_Model,'Linewidth',2,'color','k');
xlabel('$$\mathrm{\kappa}$$','interpreter','latex');
ylabel('$$\mathrm{PDF}$$','interpreter','latex');
set(gca,'TickLabelInterpreter','latex','fontsize',14);
set(gcf,'color',[1,1,1]);
leg = legend({'$\mathrm{Model}$','$\mathrm{PCE}$'});
set(leg,'Interpreter','latex');
box on;
print -depsc pdf_comp_PCE5D.eps
end




















