function [rn1,xmesh_PCE1,density_PCE1] = init_1D(eta,xpr,G,dim,L,U);

% 1D Surrogate
y = eta(:,1)'*xpr'; % active var
Y_surr = get_polyfit_surr(y,G,eta(:,1),3);

% Compute rel L-2 norm of error at MD data for kappa
np = size(xpr,1); % number of points at which MD data is available
Y_surr_data = zeros(np,1);

for i = 1:np
  Y_surr_data(i,1) = Y_surr(xpr(i,:)');
end

Nr = (sum((G-Y_surr_data).^2)).^(0.5);
Dr = (sum((G).^2)).^(0.5);
rn1 = Nr./Dr;


% SSP 1D
ssp1D(eta,xpr,G);

% Generate a million random samples in the 7D space
ts = 1e5;
samples = zeros(ts,dim);

for j = 1:dim
  samples(:,j) = unifrnd(L(j,1),U(j,1),ts,1);
end

% project those samples in [-1,1]

xrand = zeros(ts,dim);

for i = 1:ts
  for j = 1:dim
    xrand(i,j) = 2.0.*(samples(i,j)-L(j))./(U(j)-L(j)) - 1.0;
  end
end

% Estimate k using 1D surrogate for the million projected samples
Y_surr_rs = zeros(ts,1);

for i = 1:ts
  Y_surr_rs(i,1) = Y_surr(xrand(i,:)');
end

% compare pdf of kappa using 1D surrogate
pdf_comp_ssp1D(G,Y_surr_rs);

% Function Definitions

% 1D polyfit
function surr = get_polyfit_surr(y,G,eta,deg)
p = polyfit(y(:),G,deg);
surr = @(x)(polyval(p,eta'*x));
end

% plot 1D SSP
function ssp1 = ssp1D(eta,xpr,G)
figure(1)
y = eta(:,1)'*xpr';
plot(y,G,'ko','markerfacecolor','k');
xlabel('$$\mathrm{\eta_1^{T}\xi}$$','interpreter','latex','fontsize',18);
ylabel('$$\mathrm{\kappa}$$','interpreter','latex','fontsize',18);
set(gca,'TickLabelInterpreter','latex','fontsize',14);
set(gcf,'color',[1,1,1]);
box on;
print -depsc ssp1D.eps
end

% compare pdf of kappa using 1D SSP
function comp1 = pdf_comp_ssp1D(G,Y_PCE)
step = (max(Y_PCE) - min(Y_PCE)).*(1e-4);
pts_PCE = min(Y_PCE):step:max(Y_PCE);
[density_PCE1,xmesh_PCE1] = ksdensity(Y_PCE,pts_PCE);

figure;
hold on;
nbins = 15;
histogram(G,nbins,'Normalization','pdf','FaceAlpha',0.2);
plot(xmesh_PCE1,density_PCE1,'Linewidth',2,'color','b');
%plot(xmesh_Model,density_Model,'Linewidth',2,'color','k');
xlabel('$$\mathrm{\kappa}$$','interpreter','latex');
ylabel('$$\mathrm{PDF}$$','interpreter','latex');
set(gca,'TickLabelInterpreter','latex','fontsize',14);
set(gcf,'color',[1,1,1]);
leg = legend({'$\mathrm{Model}$','$\mathrm{PCE}$'});
set(leg,'Interpreter','latex');
box on;
print -depsc pdf_comp_SSP1D.eps
end

end
