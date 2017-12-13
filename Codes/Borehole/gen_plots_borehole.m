close all;
clear all;

err_samples();

function err = err_samples()
ns = [10;20;30;40;50];
err7D = [1.836;2.958e-1;1.864e-2;1.632e-2;1.934e-2];
err5D = [4.018e-1;3.136e-1;1.629e-2;2.955e-2;2.345e-3];
err4D = [1.135e-1;1.603e-1;7.151e-3;1.571e-3;1.519e-4];
err = zeros(size(ns,1),3);
err(:,1) = log10(err7D(:));
err(:,2) = log10(err5D(:));
err(:,3) = log10(err4D(:));

figure
hold on;
plot(ns,err(:,1),'-o','color','k','linewidth',2);
plot(ns,err(:,2),'-d','color','b','linewidth',2);
plot(ns,err(:,3),'-s','color','r','linewidth',2);
xlabel('$$\mathrm{Number~of~Samples}$$','interpreter','latex');
ylabel('$$\mathrm{log_{10}(\epsilon_{LOO})}$$','interpreter','latex');
xlim([8,52]);
set(gca,'fontsize',14);
set(gca,'TickLabelInterpreter','latex');
set(gcf,'color',[1,1,1]);
leg = legend({'$\mathrm{7D PCE}$','$\mathrm{5D PCE}$','$\mathrm{4D PCE}$'},'location','southwest');
set(leg,'Interpreter','latex');
box on;
print -depsc err_samples_borehole.eps
end



