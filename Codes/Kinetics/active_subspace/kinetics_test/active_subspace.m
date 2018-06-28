m = 19;
k = m + 1;
alpha = 3;
rng(123);
nsamp = floor(alpha * k * log(m));

L = zeros(1,m); U = zeros(1,m);

% Nominal values for pre-factors a1 to a19 for the 19 reactions

nom = [1.915e14,5.080e04,2.160e08,1.230e04,4.577e19,6.165e15,4.714e18,2.240e22,6.170e19,...
       6.630e13,1.690e14,1.810e13,1.450e16,3.020e12,1.202e17,1.000e13,4.820e13,9.550e06,...
       7.000e12];

L(1,:) = 0.9.*nom(1,:); U(1,:) = 1.1.*nom(1,:);

% draw samples in [-1,1]
xi = -1 + 2 * rand(nsamp, m);
%pts_xi = zeros(nsamp*(m+1),m);
%pts_xi(1:nsamp,:) = xi;
%
%% perturb the points
dxi = 2e-5;
%for i = 1:m
%  pts_xi(i*nsamp+1:(i+1)*nsamp,:) = pts_xi(1:nsamp,:);
%  pts_xi(i*nsamp+1:(i+1)*nsamp,i) = pts_xi(i*nsamp+1:(i+1)*nsamp,i) + dxi;
%end
%
%% Project points to the physical space
%pts_x = zeros(size(pts_xi,1),size(pts_xi,2));
%for i = 1:size(pts_x,1)
%  for j = 1:size(pts_x,2)
%    pts_x(i,j) = L(1,j) + 0.5*(U(1,j)-L(1,j)).*(pts_xi(i,j)+1);
%  end
%end
%
%% Save physical points to a file
%save('pts_grad.txt','pts_x','-ASCII');

id = load('record_id_grad.txt');
G = zeros(nsamp,1); Gdxi = zeros(nsamp,m);
G(:,1) = id(1:nsamp,1);
dGdxi = zeros(nsamp,m);
%
for j = 1:m
  Gdxi(:,j) = id(j*nsamp+1:(j+1)*nsamp,1);
  dGdxi(:,j) = (Gdxi(:,j) - G(:,1))./dxi;
end
%
C = zeros(m,m);
grad_f = zeros(m,1);

for i = 1:nsamp
  grad_f(:,1) = dGdxi(i,:);
  C = C + grad_f*transpose(grad_f);
end

C = (1./nsamp).*(C);
[V,D] = eig(C);

[lambda_grad, idx] = sort(diag(D), 'descend');
V = V(:,idx);

% Computing activity scores

as = zeros(m,1);

for i = 1:m
  for j=1:1
    as(i) = as(i) + lambda_grad(j).*(V(i,j).^2);
  end
end

as = as./sum(as);


figure
semilogy(abs(lambda_grad)./lambda_grad(1), '-o');
set(gca, 'fontsize', 20);
title('grad based eigs');
print -dpng eig_grad.png

eta1 = V(:,1);
eta2 = V(:,2);

%
% univariate
%
figure;
g1 = eta1'*xi';
plot(g1, G, 'ko', 'markerfacecolor', 'k')
set(gca, 'fontsize', 20);
xlabel('<eta1, x>');
ylabel('f(x)');
title('grad based SSP');
print -dpng ssp_grad.png




