rng(123);
m = 19;
%F = @(x)(borehole(x));

L = zeros(1,m); U = zeros(1,m);

% Nominal values for pre-factors a1 to a19 for the 19 reactions

nom = [1.915e14,5.080e04,2.160e08,1.230e04,4.577e19,6.165e15,4.714e18,2.240e22,6.170e19,...
       6.630e13,1.690e14,1.810e13,1.450e16,3.020e12,1.202e17,1.000e13,4.820e13,9.550e06,...
       7.000e12];

L(1,:) = 0.9.*nom(1,:); U(1,:) = 1.1.*nom(1,:);

%Algorithm 1.2 from book

k = m + 1;
alpha = 3;

N = alpha * m + 1;
M = floor( alpha *k * log(m) );

%
% get the samples
%
xi = -1 + 2 * rand(N , m); 
y = -1 + 2 * rand(M , m); 

%pts_xi = zeros(N,m);
%pts_xi(1:N,:) = xi;
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
%save('pts_grad_free.txt','pts_x','-ASCII');

id = load('record_id_grad_free.txt');
f = zeros(N,1);
f(:,1) = id(1:N);
%
% 
%%
% find the nearest p points in N for each point in M
%
p = N - 1;  %integer between m+1 and N

d_matrix = zeros(N,1);

for i=1:M
    
    for j=1:N
        
        d_matrix(j) = 0;
        
        for k=1:m
            
            d_matrix(j) = d_matrix(j) + norm(y(i,k) - xi(j,k));
               
        end
        
    end
    
    [z,index(i,:)] = sort(d_matrix);
    
    for j=1:p
        
        ip = (i-1)*p + j;
        
        points(ip,:) = xi(index(i,j),:);
        
    end
    
end

%
% formulate the least square problem
%
for np = 1 : M
    
    A = [1 points((np-1)*p+1,:)];
    
    for i = (np-1)*p+2 : np*p
        
        A = [A; 1 points(i,:)];
        
    end
    
    B = f(index(np,1));
    
    for i=2:p
        
        B = [B; f(index(np,i))];
        
    end
    
    z = A \ B;
    
    if np == 1
        
       b_matrix = z(2:m+1);
       
    else
        
       b_matrix = [b_matrix z(2:m+1)];
       
    end
    
end

%construct the covariance matrix

C = 0;

for i=1 : M
    
    z = b_matrix(:,i);
    
    C = C + z * z';
    
end

C=C/M;

[W D] = eig(C);

[lambda_loclin, idx] = sort(diag(D), 'descend');

W = W(:,idx);


eta1 = W(:,1);

eta2 = W(:,2);

% Computing activity scores

as = zeros(m,1);

for i = 1:m
  for j=1:1
    as(i) = as(i) + lambda_loclin(j).*(W(i,j).^2);
  end
end

as = as./sum(as);


figure
semilogy(lambda_loclin./lambda_loclin(1), '-o');
set(gca, 'fontsize', 20);
title('local linear eigs');
print -dpng eig_grad_free.png

%univariate SSP
figure
g1 = eta1'*xi';
plot(g1, f, 'ko', 'markerfacecolor', 'k')
set(gca, 'fontsize', 20);
xlabel('y = <eta1, x>');
ylabel('f(x)');
title('local linear SSP');
print -dpng ssp_grad_free.png
