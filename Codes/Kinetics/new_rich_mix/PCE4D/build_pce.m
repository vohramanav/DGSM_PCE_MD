close all
clearvars
uqlab

% Input data
X = dlmread(fullfile('pts_pce4D_4D_60.dat')) ;
Y = dlmread(fullfile('record_id_pce4D_19D_60.dat')) ;
MetaOpts.ExpDesign.X = X;
MetaOpts.ExpDesign.Y = Y;

L = [1.7235e14 5.5530e19 1.082e17 4.338e13];
U = [2.1065e14 6.7870e19 1.322e17 5.302e13];

nom = [1.915e14 6.170e19 1.202e17 4.820e13];

for ii = 1:4
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

%nsams = 1000;
%qoi = load('record_id_val.txt');
%E1 = verify_L2(nsams,qoi);
%save('rel_L2.txt','E1','-ASCII');
%verify_pdf(qoi,L,U);
%
%% Function Definitions
function errl = verify_L2(nsams,qoi)
pts = load('pts_val.txt');
pts4D = zeros(nsams,4);
pts4D(:,1) = pts(:,1); pts4D(:,2) = pts(:,9);
pts4D(:,3) = pts(:,15); pts4D(:,4) = pts(:,17);

Y_PCE = uq_evalModel(pts4D);
Y_Model = qoi;
NR = sum(((Y_Model-Y_PCE).^2)).^(0.5);
DR = sum((Y_Model.^2)).^(0.5);
errl = double(NR)./double(DR);
end
%
function errp = verify_pdf(qoi,L,U)
pts = zeros(1e6,4);
pts(:,1) = unifrnd(L(1),U(1),1e6,1);
pts(:,2) = unifrnd(L(2),U(2),1e6,1);
pts(:,3) = unifrnd(L(3),U(3),1e6,1);
pts(:,4) = unifrnd(L(4),U(4),1e6,1);

Y_PCE = uq_evalModel(pts);
Y_Model = qoi;
step = (max(Y_PCE) - min(Y_PCE)).*(1e-4);
pts_PCE = 0.9.*min(Y_PCE):step:1.05.*max(Y_PCE);
[density_PCE,xmesh_PCE] = ksdensity(Y_PCE,pts_PCE);

figure;
hold on;
nbins = 30;
histogram(Y_Model,nbins,'Normalization','pdf','FaceAlpha',0.2);
plot(xmesh_PCE,density_PCE,'-k','Linewidth',2);
%plot(xmesh_Model,density_Model,'Linewidth',2,'color','k');
xlabel('$$\mathrm{Ignition~Delay~(s)}$$','interpreter','latex','fontsize',20);
ylabel('$$\mathrm{PDF}$$','interpreter','latex','fontsize',20);
set(gca,'TickLabelInterpreter','latex','fontsize',18);
set(gcf,'color',[1,1,1]);
leg = legend({'$\mathrm{Model~(19D)}$','$\mathrm{PCE~(4D)}$'});
set(leg,'Interpreter','latex');
box on;
grid on;
print -depsc pdf_comp_rich.eps
end




