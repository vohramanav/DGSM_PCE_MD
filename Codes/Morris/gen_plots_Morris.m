close all
clear all

err_samples();

function Nq = nodes_order()
dim = [2;3;4];
order = linspace(1,10,10);
Nq = zeros(size(order,2),size(dim,1));
for i = 1:size(dim,1)
    Nq(:,i) = (transpose(order)+1).^(dim(i));
end

c = ['b','k','r'];
figure
hold on
for i = 1:size(dim,1)
    plot(order,Nq(:,i),'LineWidth',2,'color',c(i));
end
xlim([0,size(order,2)+1]);
xlabel('$$\mbox{Total~Order~of~Truncation}$$','interpreter','latex');
ylabel('$$\mbox{Number~of~Quadrature~Nodes}$$','interpreter','latex');
hold off
set(gca,'fontsize',14);
set(gca,'TickLabelInterpreter','latex');
set(gcf,'color',[1,1,1]);
leg = legend({'$\mathrm{d:~2}$','$\mathrm{d:~3}$','$\mathrm{d:~4}$'},'location','northwest');
set(leg,'Interpreter','latex');
box on;
print -dpng quad_comp.png
end

function err = err_samples()
ns = [5;10;15;20];
err4D = [1.2511;4.27e-1;1.27e-1;1.70e-1];
err3D = [6.54e-1;3.81e-1;8.79e-3;1.25e-20];
err2D = [3.19e-1;4.91e-26;2.49e-27;1.85e-28];
err = zeros(size(ns,1),3);
err(:,1) = log10(err4D(:));
err(:,2) = log10(err3D(:));
err(:,3) = log10(err2D(:));

figure
hold on;
plot(ns,err(:,1),'-o','color','k','linewidth',2);
plot(ns,err(:,2),'-d','color','b','linewidth',2);
plot(ns,err(:,3),'-s','color','r','linewidth',2);
xlabel('$$\mathrm{Number~of~Samples}$$','interpreter','latex');
ylabel('$$\mathrm{log_{10}(\epsilon_{LOO})}$$','interpreter','latex');
xlim([4,21]);
set(gca,'fontsize',14);
set(gca,'TickLabelInterpreter','latex');
set(gcf,'color',[1,1,1]);
leg = legend({'$\mathrm{4D PCE}$','$\mathrm{3D PCE}$','$\mathrm{2D PCE}$'},'location','southwest');
set(leg,'Interpreter','latex');
box on;
print -depsc err_samples.eps
end



