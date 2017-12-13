close all
clear all

%ub_conv();
err_samples();

function ubc = ub_conv()
data_ub = load('nsams_ub.txt');
%data_sti = load('STi_conv.mat');
col = ['c','r','b','k'];

figure
hold on;
for i = 1:size(col,2)
plot(data_ub(:,1),data_ub(:,i+1),'--o','MarkerFaceColor',col(i));
end

hold off;
xlabel('$$\mathrm{Number~of~Samples}$$','interpreter','latex');
ylabel('$$\mathrm{\hat{ub_i}}$$','interpreter','latex');
ylim([-0.01,1.1.*max(max(data_ub(:,2:end)))]);
set(gca,'fontsize',14);
set(gca,'TickLabelInterpreter','latex');
set(gcf,'color',[1,1,1]);
xlim([3,22])
leg = legend({'$\mathrm{\kappa}~[0.0764]$','$\mathrm{c}~[0.4650]$',...
              '$\mathrm{a}~[0.2298]$','$\mathrm{b}~[0.2292]$'},...
              'location','northeast');
set(leg,'Interpreter','latex');
box on;
print -depsc ub_conv_elliptic.eps
end

function err = err_samples()
filename = 'nsams_eloo_cputime.txt';
delimiterIn = ' ';
headerlinesIn = 1;
A = importdata(filename,delimiterIn,headerlinesIn);

figure
hold on;
plot(A.data(:,1),log10(A.data(:,2)),'-o','color','k','linewidth',2);
plot(A.data(:,1),log10(A.data(:,4)),'-o','color','b','linewidth',2);
ylabel('$$\mathrm{\log_{10}(\epsilon_{LOO})}$$','interpreter','latex');

yyaxis right
plot(A.data(:,1),A.data(:,3),'--o','color','k','linewidth',2);
plot(A.data(:,1),A.data(:,5),'--o','color','b','linewidth',2);
ylabel('$$\mathrm{CPU~Time}$$','interpreter','latex');

xlabel('$$\mathrm{Number~of~Training~Points}$$','interpreter','latex');
leg = legend({'$\mathrm{4D PCE}$','$\mathrm{3D PCE}$'},'location','southwest');
set(leg,'Interpreter','latex');

set(gca,'fontsize',14);
set(gca,'TickLabelInterpreter','latex');
set(gcf,'color',[1,1,1]);
box on;
grid on;
print -depsc err_samples_elliptic.eps
end
