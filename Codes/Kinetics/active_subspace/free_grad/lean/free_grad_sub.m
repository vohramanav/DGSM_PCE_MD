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
% refine xp as xpr
xpr = zeros(nsams,ncols);
xpr(:,:) = xp(1:nsams,:);

N = nsams;
p = N - 1;
d_matrix = zeros(N,1);
M = 200;
% model realizations
G = zeros(nsams,1); Gdx = zeros(nsams,dim);
G(:,1) = id_data(1:nsams,1);

% Draw M independent samples
samples = zeros(M,ncols); ypr = zeros(M,ncols);
for j = 1:dim
  samples(:,j) = unifrnd(L(j,1),U(j,1),M,1);
end

% project samples in [-1,1]
for i = 1:M
  for j = 1:dim
    ypr(i,j) = 2.0.*(samples(i,j)-L(j))./(U(j)-L(j)) - 1.0;
  end
end

for i = 1:M
  for j = 1:N
    d_matrix(j) = 0;
    for k = 1:dim
      d_matrix(j) = d_matrix(j) + norm(ypr(i,k) - xpr(j,k));
    end
  end
  [z,index(i,:)] = sort(d_matrix);
  for j = 1:p
    ip = (i-1)*p + j;
    points(ip,:) = xpr(index(i,j),:);
  end
end

% formulate the least squares problem

for np = 1:M
  A = [1 points((np-1)*p+1,:)];
  
  for i = (np-1)*p+2 : np*p        
    A = [A; 1 points(i,:)];    
  end

  B = G(index(np,1));

  for i=2:p
    B = [B; G(index(np,i))];
  end

  z = A \ B;
 
  if np == 1
    b_matrix = z(2:dim+1); %remove first entry of z 
  else
    b_matrix = [b_matrix z(2:dim+1)];
  end

end

% construct the covariance matrix

C = 0;

for i=1 : M
  z = b_matrix(:,i);
  C = C + z * z';
end

C=C/M;
[W,D] = eig(C);

[lambda, idx] = sort(diag(D), 'descend');
W = W(:,idx);

eta = zeros(dim,3);
eta(:,1) = W(:,1);
eta(:,2) = W(:,2);
eta(:,3) = W(:,3);


% Computing activity scores

as = zeros(dim,1);

for i = 1:dim
  for j=1:3
    as(i) = as(i) + lambda(j).*(W(i,j).^2);
  end
end

as = as./sum(as);

% Eigvalue plot
%eigv(lambda,dim);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[rn1,gsa_t1] = init_1D(eta,xpr,G,dim,L,U); % 1D active var analysis, 
%[rn2,xmesh_PCE2,density_PCE2,gsa_m2,gsa_t2] = init_2D(eta,xpr,G,dim,L,U,xmesh_PCE1,density_PCE1); 
% 2D active var analysis, rn: rel L2 norm of error
%rn3 = init_3D(eta,xpr,G,dim,L,U,xmesh_PCE1,density_PCE1,xmesh_PCE2,density_PCE2); % 2D active var analysis
%[rn3,gsa_m3,gsa_t3] = init_3D(eta,xpr,G,dim,L,U); % 2D active var analysis

plot_gsa(as,gsa_t1,dim);
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

function gsa = plot_gsa(as,gsa_t1,dim);
%gsa_t1 = load('gsa_t_1D.txt');
%gsa_t2 = load('gsa_t_2D.txt');
%gsa_t3 = load('gsa_t_3D.txt');
TI = zeros(dim,12);
TI(:,1) = as;
TI(:,2) = gsa_t1;
%TI(:,2) = gsa_t2;
%TI(:,3) = gsa_t3;

figure;
bar(TI);
ylabel('$$\mathrm{Activity~Scores}$$','interpreter','latex','fontsize',18);
xtickl = ({'$\mathrm{A_1}$','$\mathrm{A_2}$','$\mathrm{A_3}$','$\mathrm{A_4}$',...
           '$\mathrm{A_5}$','$\mathrm{A_6}$','$\mathrm{A_7}$','$\mathrm{A_8}$',...
           '$\mathrm{A_9}$','$\mathrm{A_{10}}$','$\mathrm{A_{11}}$','$\mathrm{A_{12}}$',...
           '$\mathrm{A_{13}}$','$\mathrm{A_{14}}$','$\mathrm{A_{15}}$','$\mathrm{A_{16}}$',...
           '$\mathrm{A_{17}}$','$\mathrm{A_{18}}$','$\mathrm{A_{19}}$','interpreter','latex'});
leg = legend({'$\mathrm{\mathcal{T}(\theta_i)}$','$\mathrm{\eta_i}$'});
set(leg,'interpreter','latex','fontsize',16,'location','Best');
set(gca,'xtick',1:dim,'xticklabel',xtickl,'fontsize',10);
set(gca,'TickLabelInterpreter','latex');
set(gcf,'color',[1,1,1]);
print -depsc as_gsa_free.eps

end

