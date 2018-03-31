function [rn2,xmesh_PCE2,density_PCE2] = init_2D(eta,xpr,G,dim,L,U,xmesh_PCE1,density_PCE1);

% 2D Polynomial surface fit
y1 = eta(:,1)'*xpr';
y2 = eta(:,2)'*xpr';
sf = fit([y1',y2'],G,'poly55');

% compute rel L-2 norm of error
np = size(xpr,1);
Y_surr_data = zeros(np,1);

for i = 1:np
  Y_surr_data(i) = sf(y1(i),y2(i));
end

Nr = (sum((G-Y_surr_data).^2)).^(0.5);
Dr = (sum((G).^2)).^(0.5);
rn2 = Nr./Dr; 

% SSP 2D
%ssp2D(y1,y2,G,sf); % 3D plot with surface fit
ssp2D_2(y1,y2,G); % 2D plot showing trends in k along first two active variables

% compare pdf with the histogram'
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

Y_surr_rs = zeros(ts,1);
y1_rs = eta(:,1)'*xrand';
y2_rs = eta(:,2)'*xrand';

for i = 1:ts
  Y_surr_rs(i) = sf(y1_rs(i),y2_rs(i));
end

%pdf_comp_ssp2D(G,Y_surr_rs);

% Function Definitions

function ssp2 = ssp2D(y1,y2,G,sf)
figure;
%hold on;
%scatter3(y1,y2,G,'ko','markerfacecolor','k');
plot(sf,[y1',y2'],G);
shading interp;
xlabel('$$\mathrm{\eta_1^{T}\xi}$$','interpreter','latex','fontsize',18);
ylabel('$$\mathrm{\eta_2^{T}\xi}$$','interpreter','latex','fontsize',18);
zlabel('$$\mathrm{\kappa}$$','interpreter','latex','fontsize',18);
set(gca,'TickLabelInterpreter','latex','fontsize',14);
set(gcf,'color',[1,1,1]);
box on;
print -depsc ssp2D.eps
end

function ssp2_2 = ssp2D_2(y1,y2,G);
figure;
ps = 100; % size of the point
scatter(y1,y2,ps,G,'filled');
shading interp;
xlabel('$$\mathrm{\eta_1^{T}\xi}$$','interpreter','latex','fontsize',18);
ylabel('$$\mathrm{\eta_2^{T}\xi}$$','interpreter','latex','fontsize',18);
set(gca,'TickLabelInterpreter','latex','fontsize',14);
set(gcf,'color',[1,1,1]);
box on;
c = colorbar();
set(c,'TickLabelInterpreter','latex');
print -depsc ssp2D_2.eps
end


% compare pdf of kappa using 2D SSP
function comp2 = pdf_comp_ssp2D(G,Y_PCE)
step = (max(Y_PCE) - min(Y_PCE)).*(1e-4);
pts_PCE = min(Y_PCE):step:max(Y_PCE);
[density_PCE2,xmesh_PCE2] = ksdensity(Y_PCE,pts_PCE);

figure;
hold on;
nbins = 15;
histogram(G,nbins,'Normalization','pdf','FaceAlpha',0.2);
plot(xmesh_PCE1,density_PCE1,'Linewidth',1.5,'color','k');
plot(xmesh_PCE2,density_PCE2,'Linewidth',1.5,'color','b');
%plot(xmesh_Model,density_Model,'Linewidth',2,'color','k');
xlabel('$$\mathrm{\kappa}$$','interpreter','latex');
xlim([0,50]);
ylabel('$$\mathrm{PDF}$$','interpreter','latex');
set(gca,'TickLabelInterpreter','latex','fontsize',14);
set(gcf,'color',[1,1,1]);
%leg = legend({'$\mathrm{MD}$','$\mathrm{2D~Subspace}$'});
leg = legend({'$\mathrm{MD}$','$\mathrm{1D~Subspace}$','$\mathrm{2D~Subspace}$'});
set(leg,'Interpreter','latex');
box on;
print -depsc pdf_comp_SSP2D.eps
end


end
