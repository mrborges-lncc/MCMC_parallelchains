clear;
close all
loc=300;
jump=8;
N=00;
B=250.;
A=0;
% 1 m3 = 6.2898105697751 bbl
fat = 6.2898105697751;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Nch_ini = 0;
Nch_fim = 0;
Nchains = Nch_fim - Nch_ini + 1;
Nini = repmat(0, 1, Nchains);
Nfim = [1999];
Nfim = Nfim(Nch_ini+1:Nch_fim+1);
Nt   = (Nfim-Nini)+1;
chains = [Nch_ini:1:Nch_fim];
base_name = 'prod_TwoPhase3D_KC_MC'
% base_name = 'prod_TwoPhase3DMC_KC'
hom  = '~/Dropbox/PROJETO_MCMC_RIGID/MCMC_parallelchains/';
homf = '~/Dropbox/PROJETO_MCMC_RIGID/paper/figuras/';
dados=load([hom 'MonteCarlo/twophaseflow/exp/prod/prod_referencia_0.dat']);
%dados=load(['~/MCMC_parallelchains/twophaseflow/exp/prod/prod_referencia_0.dat']);
ref=dados;
home = [hom 'MonteCarlo/twophaseflow/exp/prod/'];
%
data=zeros(size(ref,1),size(ref,2),sum(Nt));
k = 0;
for j=1:Nchains
    n = num2str(chains(j),'%d');
    for i=Nini(j):Nfim(j)
        k = k+1;
        istr=num2str(i,5);
        fprintf('\nSample no. %s',istr);
        file_name   = [home base_name '_' istr '.dat'];
        data(:,:,k) = load(file_name);
%         close all
%         figure(33)
%         plot(data(:,1,k),data(:,2,k),'-k','LineWidth',3);
%         hold on
%         plot(data(:,1,k),data(:,3,k),'-r','LineWidth',3);
%         plot(data(:,1,k),data(:,4,k),'-b','LineWidth',3);
%         plot(data(:,1,k),data(:,5,k),'-g','LineWidth',3);
%         ylim([0 300]);
    end
end
dmedio = mean(data,3);
erro   = std(data,0,3);
%
% Create figure
figure1 = figure()

C=min(dados(:,1));
D=max(dados(:,1))*1.01;

dasp=[1 1.*(B-A)/(D-C) 200];
% Create axes
axes1 = axes('Parent',figure1,'FontSize',14,'FontName','Times New Roman',...
    'DataAspectRatio',dasp);
box(axes1,'on');
hold(axes1,'all');

dados = dmedio;
errorbar(dados(1:jump:end,1),dados(1:jump:end,2),erro(1:jump:end,2),...
    'Parent',axes1,'Color',[0.85 0.33 0.10],'MarkerSize',6,'Marker','o',...
    'LineStyle','none','DisplayName','mean 1','LineWidth',0.5)
errorbar(dados(1:jump:end,1),dados(1:jump:end,3),erro(1:jump:end,3),...
    'Parent',axes1,'Color',[0.07 0.62 1],'MarkerSize',6,'Marker','^',...
    'LineStyle','none','DisplayName','mean 2','LineWidth',0.5)
errorbar(dados(1:jump:end,1),dados(1:jump:end,4),erro(1:jump:end,4),...
    'Parent',axes1,'Color',[0.93 0.69 0.13],'MarkerSize',6,'Marker','s',...
    'LineStyle','none','DisplayName','mean 3','LineWidth',0.5)
errorbar(dados(1:jump:end,1),dados(1:jump:end,5),erro(1:jump:end,4),...
    'Parent',axes1,'Color',[0 0 0],'MarkerSize',6,'Marker','v',...
    'LineStyle','none','DisplayName','mean 4','LineWidth',0.5)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dados=ref;
plot(dados(2:end,1),dados(2:end,2),'Parent',axes1,'Color',[0.85 0.33 0.10],...
    'MarkerSize',6,'LineWidth',2,'DisplayName','ref.  1')
plot(dados(2:end,1),dados(2:end,3),'Parent',axes1,'Color',[0.07 0.62 1],...
    'MarkerSize',6,'LineWidth',2,'DisplayName','ref.  2')
plot(dados(2:end,1),dados(2:end,4),'Parent',axes1,'Color',[0.93 0.69 0.13],...
    'MarkerSize',6,'LineWidth',2,'DisplayName','ref.  3')
plot(dados(2:end,1),dados(2:end,5),'Parent',axes1,'Color',[0 0 0],...
    'MarkerSize',6,'LineWidth',2,'DisplayName','ref.  4')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plot([loc loc],[A B],'Parent',axes1,'Color',[0 0 0],...
    'MarkerSize',6,'LineWidth',1,'LineStyle','--',...
    'DisplayName','$t_s$')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xlim(axes1,[0 D])
ylim(axes1,[A B])
% Create xlabel
xlabel('$t (day)$','FontSize',16,'FontName','Times New Roman',...
    'FontAngle','italic','Interpreter','latex');

% Create ylabel
ylabel('Oil rate ($m^3/day$)','FontSize',16,'FontName',...
    'Times New Roman','FontAngle','italic','Interpreter','latex');

% Set the remaining axes properties
set(axes1,'FontName',...
    'Times New Roman','FontSize',14,'TickDir','both','TickLabelInterpreter',...
    'latex','XMinorTick','on','YMinorTick','on');

% Create legend
legend1 = legend(axes1,'show');
set(legend1,'Location','NorthEast','FontSize',11,'Interpreter','latex');
set(legend1,'Box','off');

% Print
base=[homf base_name];
%print('-djpeg90',base)
print('-depsc','-r600',base)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Production
sz = size(data);
producao = zeros(sz(1),2,sz(3));
for k = 1:size(data,3)
    t = data(:,1,k);
    p = sum(data(:,2:end,k),2);
    pa= 0.0;
    prod = [0.0 0.0];
    for j = 2:size(t,1)
        dt = t(j) - t(j-1);
        pr = 0.5*(p(j) + p(j-1));
        pa = pa + pr * dt;
        prod = [prod; t(j) pa];
    end
    producao(:,:,k) = prod;
end
t = ref(:,1);
p = sum(ref(:,2:end),2);
prodref = [0.0 0.0];
pa= 0.0;
for j = 2:size(t,1)
    dt = t(j) - t(j-1);
    pr = 0.5*(p(j) + p(j-1));
    pa = pa + pr * dt;
    prodref = [prodref; t(j) pa];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create figure
clf; close all
figure2 = figure(2)

C=min(min(prodref(:,1,:)));
D=max(max(prodref(:,1,:)))*1.01;
A=min(min(prodref(:,2,:)));
B=max(max(prodref(:,2,:)))*1.2;

dasp=[1 1.*(B-A)/(D-C) 200];
% Create axes
axes1 = axes('Parent',figure2,'FontSize',14,'FontName','Times New Roman',...
    'DataAspectRatio',dasp);
box(axes1,'on');
hold(axes1,'all');
for k = 1:size(data,3)
    h=plot(producao(1:end,1,k),producao(1:end,2,k),'Parent',axes1,...
        'Color',[0.65 0.65 0.65],'MarkerSize',6,'LineWidth',0.5,...
        'DisplayName','');
    h.Annotation.LegendInformation.IconDisplayStyle = 'off';
end
plot(prodref(1:end,1),prodref(1:end,2),'Parent',axes1,...
    'Color',[1 0 0],'MarkerSize',6,'LineWidth',3,...
    'DisplayName','ref.');
prod = mean(producao,3);
plot(prod(1:end,1),prod(1:end,2),'Parent',axes1,...
    'Color',[0 0 0],'MarkerSize',6,'LineWidth',2,...
    'DisplayName','mean');
sprod = std(producao,0,3);
errorbar(prod(1:jump:end,1),prod(1:jump:end,2),sprod(1:jump:end,2),...
    'Parent',axes1,'Color',[0. 0. 0.],'MarkerSize',4,'Marker','none',...
    'LineStyle','-','DisplayName','error','LineWidth',0.5);
xlim(axes1,[C D]);
ylim(axes1,[A B]);
% Create xlabel
xlabel('$t (day)$','FontSize',16,'FontName','Times New Roman',...
    'FontAngle','italic','Interpreter','latex');

% Create ylabel
ylabel('Cumulated oil production ($m^3$)','FontSize',16,'FontName',...
    'Times New Roman','FontAngle','italic','Interpreter','latex');

% Set the remaining axes properties
set(axes1,'FontName',...
    'Times New Roman','FontSize',14,'TickDir','both','TickLabelInterpreter',...
    'latex','XMinorTick','on','YMinorTick','on');

% Create legend
legend1 = legend(axes1,'show');
set(legend1,'Location','NorthWest','FontSize',10);
set(legend1,'Box','off');

% Print
base=[homf 'total_' base_name];
%print('-djpeg90',base)
print('-depsc','-r300',base)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
