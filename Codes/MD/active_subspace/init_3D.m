function rn3 = init_3D(eta,xpr,G,dim,L,U,xmesh_PCE1,density_PCE1,xmesh_PCE2,density_PCE2);

% 3D Polynomial surface fit

y1 = eta(:,1)'*xpr';
y2 = eta(:,2)'*xpr';
y3 = eta(:,3)'*xpr';

s = size(y1,2);
Y = zeros(s,3);
Y(:,1) = y1;
Y(:,2) = y2;
Y(:,3) = y3;

reg=MultiPolyRegress(Y,G,3);

% compute rel L-2 norm of error
np = size(xpr,1);
Y_surr_data = zeros(np,1);

for i = 1:np
  np=Y(i,:);
  ns=repmat(np,[length(reg.PowerMatrix) 1]).^reg.PowerMatrix;
  es=ones(length(reg.PowerMatrix),1);
  for ii=1:size(reg.PowerMatrix,2)
    es=es.*ns(:,ii);
  end
  Y_surr_data(i) = reg.Coefficients'*es; % The estimate for the new data point.
end

Nr = (sum((G-Y_surr_data).^2)).^(0.5);
Dr = (sum((G).^2)).^(0.5);
rn3 = Nr./Dr;

% compare pdf with histogram

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
y3_rs = eta(:,3)'*xrand';
sn = size(y1_rs,2);
Y_rs = zeros(sn,3);
Y_rs(:,1) = y1_rs;
Y_rs(:,2) = y2_rs;
Y_rs(:,3) = y3_rs;

for i = 1:ts
  np=Y_rs(i,:);
  %np = [-0.4716   -0.2264    0.6605];
  ns=repmat(np,[length(reg.PowerMatrix) 1]).^reg.PowerMatrix;
  es=ones(length(reg.PowerMatrix),1);
  for ii=1:size(reg.PowerMatrix,2)
    es=es.*ns(:,ii);
  end
  Y_surr_rs(i) = reg.Coefficients'*es; % The estimate for the new data point.
end

pdf_comp_ssp3D(G,Y_surr_rs);

% compare pdf of kappa using 2D SSP
function comp3 = pdf_comp_ssp3D(G,Y_PCE)
step = (max(Y_PCE) - min(Y_PCE)).*(1e-4);
pts_PCE = min(Y_PCE):step:max(Y_PCE);
[density_PCE3,xmesh_PCE3] = ksdensity(Y_PCE,pts_PCE);

figure;
hold on;
nbins = 15;
histogram(G,nbins,'Normalization','pdf','FaceAlpha',0.2);
plot(xmesh_PCE1,density_PCE1,'Linewidth',1.5,'color','k');
plot(xmesh_PCE2,density_PCE2,'Linewidth',1.5,'color','b');
plot(xmesh_PCE3,density_PCE3,'Linewidth',1.5,'color','r');
%plot(xmesh_Model,density_Model,'Linewidth',2,'color','k');
xlabel('$$\mathrm{\kappa}$$','interpreter','latex');
xlim([0,50]);
ylabel('$$\mathrm{PDF}$$','interpreter','latex');
set(gca,'TickLabelInterpreter','latex','fontsize',14);
set(gcf,'color',[1,1,1]);
%leg = legend({'$\mathrm{MD}$','$\mathrm{2D~Subspace}$'});
leg = legend({'$\mathrm{MD}$','$\mathrm{1D~Subspace}$',...
              '$\mathrm{2D~Subspace}$','$\mathrm{3D~Subspace}$'});
set(leg,'Interpreter','latex');
box on;
print -depsc pdf_comp_SSP3D.eps
end





end
