function [lambda, W, x, f, y] = compute_active_subspace(F, m, dist);


close all;

k = m + 1;
alpha = 5;

nsamples = floor(alpha * k * log(m));


%
% get the samples
%
switch dist
   case 'uniform'
      x = -1 + 2 * rand(nsamples , m);

   case 'gaussian'
      x = randn(nsamples , m);

   otherwise
       error('unknown option, use uniform or gaussian');
end

fprintf('running %i samples\n', nsamples);
   
parfor i = 1 : nsamples
   fprintf('%i\n', i);
   [f(i) df(i,:)] = F(x(i,:)');
end

C = 0;

for i = 1 : nsamples
   dfi = df(i,:)';
   C = C + dfi * dfi'; 
end
C = C / nsamples;

[W D] = eig(C);
[lambda, idx] = sort(diag(D), 'descend');
W = W(:,idx);


eta1 = W(:,1);
eta2 = W(:,2);

%
% plot of eigenvalues
%
figure();
semilogy(1:8,lambda./lambda(1), '-o','MarkerFaceColor','k');
xlabel('$$\mathrm{\mbox{Index}}$$','interpreter','latex','fontsize',20);
ylabel('$$\mathrm{\mbox{Eigenvalue}}$$','interpreter','latex','fontsize',20);
set(gca,'ticklabelinterpreter','latex','fontsize',18);
box on;
grid on;
print -depsc lam.eps


%
% univariate SSP
%
figure;
y = eta1'*x';
plot(y, f, 'ko', 'markerfacecolor', 'k')
xlabel('$$\mathrm{\eta_1^\top x}$$','interpreter','latex','fontsize',20);
ylabel('$$\mathrm{I(t^\ast)}$$','interpreter','latex','fontsize',20);
set(gca,'ticklabelinterpreter','latex','fontsize',18);
box on;
print -depsc uni_SSP.eps

%
% bivariate SSP
%
%figure
%g1 = y;
%g2 = eta2'*x';
%scatter3(g1, g2, f, 'ko', 'markerfacecolor', 'k');
%set(gca, 'fontsize', 20);

return 




%
% surrogate model stuff
%
p = polyfit(y(:),f(:),2)

W1 = W(:,1);
F = @(x)(polyval(p,W1'*x));

