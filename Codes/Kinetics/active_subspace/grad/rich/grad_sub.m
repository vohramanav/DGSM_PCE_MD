close all;
clear all;
rng(100);
set(0,'DefaultFigureVisible','off');

dim = 19;
dX = zeros(dim,1); 
L = zeros(dim,1); U = zeros(dim,1);
id_data = load('record_id_as.txt');
params = load('pts_as.txt');

% Nominal values of the parameters
nom = [1.915e14,5.080e04,2.160e08,1.230e04,4.577e19,6.165e15,4.714e18,2.240e22,6.170e19,...
       6.630e13,1.690e14,1.810e13,1.450e16,3.020e12,1.202e17,1.000e13,4.820e13,9.550e06,...
       7.000e12];

L(:,1) = 0.9.*nom(:,1); U(:,1) = 1.1.*nom(:,1);

% Project params in [-1,1]

nrows = size(params,1); % number of rows
ncols = size(params,2); % number of cols

xp = zeros(nrows,ncols);

for i = 1:nrows
  for j = 1:ncols
    xp(i,j) = 2.0.*(params(i,j)-L(j))./(U(j)-L(j)) - 1.0;
  end
end

nsams = 60;


for j = 1:dim
  dX(j) = 1e-3.*(U(j,1)-L(j,1));
end

G = zeros(nsams,1); Gdx = zeros(nsams,dim);
G(:,1) = id_data(1:nsams,1);
dGdx = zeros(nsams,dim);

for j = 1:dim
  Gdx(:,j) = id_data(j*nsams+1:(j+1)*nsams,1);
  dGdx(:,j) = (Gdx(:,j) - G(:,1))./dX(j);
end

% Computing the C matrix (Eq. 1.2)

C = zeros(dim,dim);
grad_f = zeros(dim,1);

for i = 1:nsams
  grad_f(:,1) = dGdx(i,:);
  C = C + grad_f*transpose(grad_f);
end

C = (1./nsams).*(C);
[W,D] = eig(C);

[lambda, idx] = sort(diag(D), 'descend');
W = W(:,idx);

eta = zeros(dim,3);
eta(:,1) = W(:,1);
eta(:,2) = W(:,2);
eta(:,3) = W(:,3);
data = load('W_free.txt');

figure;
hold on
plot(W(:,1),'o');
plot(data(:,1),'*');
print -dpng compare.png


% Computing activity scores

as = zeros(dim,1);

for i = 1:dim
  for j=1:1
    as(i) = as(i) + lambda(j).*(W(i,j).^2);
  end
end

as = as./sum(as);

% Eigvalue plot
%eigv(lambda,dim);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%[rn1,y] = init_1D(eta,xp(1:nsams,:),G,dim,L,U); % 1D active var analysis, 
%[rn2,xmesh_PCE2,density_PCE2,gsa_m2,gsa_t2] = init_2D(eta,xpr,G,dim,L,U,xmesh_PCE1,density_PCE1); 
% 2D active var analysis, rn: rel L2 norm of error
%rn3 = init_3D(eta,xpr,G,dim,L,U,xmesh_PCE1,density_PCE1,xmesh_PCE2,density_PCE2); % 2D active var analysis
%[rn3,gsa_m3,gsa_t3] = init_3D(eta,xpr,G,dim,L,U); % 2D active var analysis

%plot_gsa(as,dim);
%plot_rn(rn1,rn2,rn3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function eigen = eigv(lambda,dim)
figure(1)
%plot(1:7,e,'--o','MarkerFaceColor','b');
semilogy(lambda./lambda(1),'-o');
xlabel('$$\mathrm{index}$$','interpreter','latex','fontsize',18);
ylabel('$$\mathrm{Eigenvalues}$$','interpreter','latex','fontsize',18);
set(gca, 'xtick',1:dim,'fontsize',14);
set(gca,'TickLabelInterpreter','latex');
set(gcf,'color',[1,1,1]);
box on;
grid on;
print -depsc eigv.eps
end

function prn = plot_rn(rn1,rn2,rn3);
figure;
plot([1 2 3],[rn1 rn2 rn3],'o--','MarkerFaceColor','k');
xlabel('$$\mathrm{Active~Subspace~Dimension}$$','interpreter','latex','fontsize',18);
ylabel('$$\mathrm{Relative~L-2~Norm~of~Error}$$','interpreter','latex','fontsize',18);
xticks([1 2 3]);
xlim([0,4]);
set(gca,'TickLabelInterpreter','latex');
set(gcf,'color',[1,1,1]);
box on;
print -depsc rn.eps
end

function gsa = plot_gsa(as,dim);
%gsa_t1 = load('gsa_t_1D.txt');
%gsa_t2 = load('gsa_t_2D.txt');
%gsa_t3 = load('gsa_t_3D.txt');
TI = zeros(dim,1);
TI(:,1) = as;
%TI(:,2) = gsa_t1;
%TI(:,2) = gsa_t2;
%TI(:,3) = gsa_t3;

figure;
bar(TI,'BarWidth',0.5);
ylabel('$$\mathrm{Activity~Scores}$$','interpreter','latex','fontsize',18);
xtickl = ({'$\mathrm{A_1}$','$\mathrm{A_2}$','$\mathrm{A_3}$','$\mathrm{A_4}$',...
           '$\mathrm{A_5}$','$\mathrm{A_6}$','$\mathrm{A_7}$','$\mathrm{A_8}$',...
           '$\mathrm{A_9}$','$\mathrm{A_{10}}$','$\mathrm{A_{11}}$','$\mathrm{A_{12}}$',...
           '$\mathrm{A_{13}}$','$\mathrm{A_{14}}$','$\mathrm{A_{15}}$','$\mathrm{A_{16}}$',...
           '$\mathrm{A_{17}}$','$\mathrm{A_{18}}$','$\mathrm{A_{19}}$','interpreter','latex'});
%leg = legend({'$\mathrm{\mathcal{T}(\theta_i)}$','$\mathrm{\eta_i}$'});
%set(leg,'interpreter','latex','fontsize',16,'location','Best');
set(gca,'xtick',1:dim,'xticklabel',xtickl,'fontsize',10);
%set(leg,'interpreter','latex','fontsize',14);
set(gca,'TickLabelInterpreter','latex');
set(gcf,'color',[1,1,1]);
print -depsc as.eps

figure;
bar(as,'BarWidth',0.5);
ylabel('$$\mathrm{Activity~Scores}$$','interpreter','latex','fontsize',18);
xtickl = ({'$$\mathrm{A}$$','$$\mathrm{B}$$','$$\mathrm{p}$$','$$\mathrm{q}$$',...
           '$$\mathrm{\alpha}$$','$$\mathrm{\lambda}$$','$$\mathrm{\gamma}$$',...
           'interpreter','latex'});
set(gca,'xtick',1:7,'xticklabel',xtickl,'fontsize',16);
set(gca,'TickLabelInterpreter','latex');
set(gcf,'color',[1,1,1]);
print -depsc activity_scores.eps

end

