close all
clear all

err_samples();

function err = err_samples()
filename = 'nsams_eloo.txt';
delimiterIn = ' ';
headerlinesIn = 1;
A = importdata(filename,delimiterIn,headerlinesIn);

figure
semilogy(A.data(:,1),A.data(:,2),'-o','color','k','linewidth',2,'MarkerFaceColor','k');
hold on;
semilogy(A.data(:,1),A.data(:,3),'--o','color','b','linewidth',2,'MarkerFaceColor','b');
ylabel('$$\mathrm{\log_{10}(\epsilon_{LOO})}$$','interpreter','latex','fontsize',20);

%yyaxis right
%plot(A.data(:,1),A.data(:,3),'--o','color','k','linewidth',2);
%plot(A.data(:,1),A.data(:,5),'--o','color','b','linewidth',2);
%ylabel('$$\mathrm{CPU~Time}$$','interpreter','latex');

xlabel('$$\mathrm{Number~of~Training~Points}$$','interpreter','latex','fontsize',20);
xlim([0,70]);
leg = legend({'$\mathrm{19D~PCE}$','$\mathrm{4D~PCE}$'},'location','southwest');
set(leg,'Interpreter','latex');

set(gca,'fontsize',18);
set(gca,'TickLabelInterpreter','latex');
set(gcf,'color',[1,1,1]);
box on;
grid on;
print -depsc err_samples_kinetics.eps
end
