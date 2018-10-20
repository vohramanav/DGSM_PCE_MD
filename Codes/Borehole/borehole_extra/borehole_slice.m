function [v1 v2] = borehole_slice(i)

x = [-1 : 0.01 : 1];
nx = length(x);
xi = zeros(7,1);
for j = 1 : nx
   xi(i) = x(j); 
   y(j) = borehole(xi);
end
close all
figure(1); 
plot(x, y, 'linewidth',2);

figure
dy = gradient(y);
semilogy(x, dy.^2, 'linewidth',2);

figure;
plot(x, dy, 'linewidth',2);



x = -1 + 2 * rand(1e4,1);
nx = length(x);
for j = 1 : nx
   xi = zeros(7,1); xi(i) = x(j);
   [f Df] = borehole(xi);
   yy(j) = f;
   dyy(j) = abs(Df(i));
end
v1 = std(yy);
v2 = std(dyy);
disp([v1 v2])
figure
bar([v1 v2])
set(gca, 'yscale', 'log')

figure; 
subplot(121);
hist(yy);
subplot(122);
hist(dyy);
