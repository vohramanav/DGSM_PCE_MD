function [rn1] = init_1D(eta,xpr,G,dim,L,U);

% 1D Surrogate
y = eta(:,1)'*xpr'; % active var
Y_surr = get_polyfit_surr(y,G,eta(:,1),2);

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
%ts = 1e5;
%samples = zeros(ts,dim);
%
%for j = 1:dim
%  samples(:,j) = unifrnd(L(j,1),U(j,1),ts,1);
%end
%
%% project those samples in [-1,1]
%
%xrand = zeros(ts,dim);
%
%for i = 1:ts
%  for j = 1:dim
%    xrand(i,j) = 2.0.*(samples(i,j)-L(j))./(U(j)-L(j)) - 1.0;
%  end
%end
%
%% Estimate k using 1D surrogate for the million projected samples
%Y_surr_rs = zeros(ts,1);
%
%for i = 1:ts
%  Y_surr_rs(i,1) = Y_surr(xrand(i,:)');
%end
%
%% compare pdf of kappa using 1D surrogate
%pdf_comp_ssp1D(G,Y_surr_rs);
%[gsa_t1] = sobol(L,U,dim);

% Function Definitions

% 1D polyfit
function surr = get_polyfit_surr(y,G,eta,deg)
p = polyfit(y(:),G,deg);
surr = @(x)(polyval(p,eta'*x));
end

function [gsa_t1] = sobol(L,U,dim)
N = 1e5; % number of samples
A = zeros(N,dim); B = zeros(N,dim);

for i = 1:dim
  A(:,i) = unifrnd(L(i,1),U(i,1),N,1);
  B(:,i) = unifrnd(L(i,1),U(i,1),N,1);
end

% project those samples in [-1,1]

An = zeros(N,dim); Bn = zeros(N,dim);

for j = 1:N
  for i = 1:dim
    An(j,i) = 2.0.*(A(j,i)-L(i))./(U(i)-L(i)) - 1.0;
    Bn(j,i) = 2.0.*(B(j,i)-L(i))./(U(i)-L(i)) - 1.0;
  end
end

% constructing the modified matrices A(B)i and B(A)i
An_Bn_1 = An; An_Bn_1(:,1) = Bn(:,1);
An_Bn_2 = An; An_Bn_2(:,2) = Bn(:,2);
An_Bn_3 = An; An_Bn_3(:,3) = Bn(:,3);
An_Bn_4 = An; An_Bn_4(:,4) = Bn(:,4);
An_Bn_5 = An; An_Bn_5(:,5) = Bn(:,5);
An_Bn_6 = An; An_Bn_6(:,6) = Bn(:,6);
An_Bn_7 = An; An_Bn_7(:,7) = Bn(:,7);
An_Bn_8 = An; An_Bn_8(:,8) = Bn(:,8);
An_Bn_9 = An; An_Bn_9(:,9) = Bn(:,9);
An_Bn_10 = An; An_Bn_10(:,10) = Bn(:,10);
An_Bn_11 = An; An_Bn_11(:,11) = Bn(:,11);
An_Bn_12 = An; An_Bn_12(:,12) = Bn(:,12);
An_Bn_13 = An; An_Bn_13(:,13) = Bn(:,13);
An_Bn_14 = An; An_Bn_14(:,14) = Bn(:,14);
An_Bn_15 = An; An_Bn_15(:,15) = Bn(:,15);
An_Bn_16 = An; An_Bn_16(:,16) = Bn(:,16);
An_Bn_17 = An; An_Bn_17(:,17) = Bn(:,17);
An_Bn_18 = An; An_Bn_18(:,18) = Bn(:,18);
An_Bn_19 = An; An_Bn_19(:,19) = Bn(:,19);

%surrogate evaluations
fAn = zeros(N,1); 
fAn_Bn = zeros(N,dim);

for i = 1:N
  fAn(i,1) = Y_surr(An(i,:)');
  fAn_Bn(i,1) = Y_surr(An_Bn_1(i,:)');  
  fAn_Bn(i,2) = Y_surr(An_Bn_2(i,:)');  
  fAn_Bn(i,3) = Y_surr(An_Bn_3(i,:)'); 
  fAn_Bn(i,4) = Y_surr(An_Bn_4(i,:)'); 
  fAn_Bn(i,5) = Y_surr(An_Bn_5(i,:)'); 
  fAn_Bn(i,6) = Y_surr(An_Bn_6(i,:)');  
  fAn_Bn(i,7) = Y_surr(An_Bn_7(i,:)'); 
  fAn_Bn(i,8) = Y_surr(An_Bn_8(i,:)'); 
  fAn_Bn(i,9) = Y_surr(An_Bn_9(i,:)'); 
  fAn_Bn(i,10) = Y_surr(An_Bn_10(i,:)'); 
  fAn_Bn(i,11) = Y_surr(An_Bn_11(i,:)'); 
  fAn_Bn(i,12) = Y_surr(An_Bn_12(i,:)'); 
  fAn_Bn(i,13) = Y_surr(An_Bn_13(i,:)'); 
  fAn_Bn(i,14) = Y_surr(An_Bn_14(i,:)'); 
  fAn_Bn(i,15) = Y_surr(An_Bn_15(i,:)'); 
  fAn_Bn(i,16) = Y_surr(An_Bn_16(i,:)'); 
  fAn_Bn(i,17) = Y_surr(An_Bn_17(i,:)'); 
  fAn_Bn(i,18) = Y_surr(An_Bn_18(i,:)'); 
  fAn_Bn(i,19) = Y_surr(An_Bn_19(i,:)'); 
end

f_total = zeros((dim+1).*N,1);
f_total(1:N,1) = fAn;
f_total(N+1:end,1) = reshape(fAn_Bn,N*dim,1);
f0_t = mean(f_total); Dr_t = var(f_total);

Nr_t = zeros(dim,1); gsa_t = zeros(dim,1);

for i = 1:dim
   diff = fAn - fAn_Bn(:,i);
   Nr_t(i) = (1.0/N).*dot(fAn,diff); 
end

gsa_t1 = Nr_t./Dr_t;

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
leg = legend({'$\mathrm{MD~Prediction}$','$\mathrm{1D~Subspace}$'});
set(leg,'Interpreter','latex');
box on;
print -depsc free_pdf_comp_SSP1D.eps
end

end
