clear;
close all
loc=100;
jump=4;
N=60;
B=230;
A=70;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Nch_ini = 0;
Nch_fim = 0;
Nchains = Nch_fim - Nch_ini + 1;
Nini = repmat(0, 1, Nchains);
Nfim = [1999];
Nfim = Nfim(Nch_ini+1:Nch_fim+1);
Nt   = (Nfim-Nini)+1;
chains = [Nch_ini:1:Nch_fim];
base_name = 'presinj_TwoPhase3DMC_only_perm'
base_name = 'presinj_TwoPhase3DMC_KC'
hom  = '~/Dropbox/PROJETO_MCMC_RIGID/MCMC_parallelchains/';
homf = '~/Dropbox/PROJETO_MCMC_RIGID/paper/figuras/';
dados=load([hom 'MonteCarlo/twophaseflow/exp000/pres/pres_referencia_0.dat']);
ref=dados;
home = [hom 'MonteCarlo/twophaseflow/exp000/pres/'];
%
data=zeros(size(ref,1),size(ref,2),sum(Nt));
k = 0;
for j=1:Nchains
    n = num2str(chains(j),'%d');
    for i=Nini(j):Nfim(j)
        k = k+1;
        istr=num2str(i,5);
        file_name   = [home base_name '_' istr '.dat']
        data(:,:,k) = load(file_name);
    end
end
dmedio = mean(data,3);
erro   = std(data,0,3);
%
% Create figure
figure1 = figure()
%B=1e-2;
%A=1e-5;
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
    'Parent',axes1,'Color',[1 0 0],'MarkerSize',4,'Marker','o',...
    'LineStyle','none','DisplayName','mean','LineWidth',0.5)
% errorbar(dados(1:jump:end,1),dados(1:jump:end,3),erro(1:jump:end,3),...
%     'Parent',axes1,'Color',[0 0 1],'MarkerSize',4,'Marker','s',...
%     'LineStyle','none','DisplayName','mean 2','LineWidth',1)
% errorbar(dados(1:jump:end,1),dados(1:jump:end,4),erro(1:jump:end,4),...
%     'Parent',axes1,'Color',[0 0 0],'MarkerSize',4,'Marker','o',...
%     'LineStyle','none','DisplayName','mean 3','LineWidth',1)
% errorbar(dados(1:jump:end,1),dados(1:jump:end,5),erro(1:jump:end,4),...
%     'Parent',axes1,'Color',[0 1 0],'MarkerSize',4,'Marker','s',...
%     'LineStyle','none','DisplayName','mean 4','LineWidth',1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dados=ref;
%dados=load('../../MCMC/reject_prod/prod_chn0-20_16000.dat');
%dados=load('../conc/conc_amostra_0.dat');
plot(dados(2:end,1),dados(2:end,2),'Parent',axes1,'Color',[1 0 0],...
    'MarkerSize',6,'LineWidth',2,'DisplayName','ref.')
% plot(dados(:,1),dados(:,3),'Parent',axes1,'Color',[0 0 1],...
%     'MarkerSize',6,'LineWidth',2,'DisplayName','ref.  2')
% plot(dados(:,1),dados(:,4),'Parent',axes1,'Color',[0 0 0],...
%     'MarkerSize',6,'LineWidth',2,'DisplayName','ref.  3')
% plot(dados(:,1),dados(:,5),'Parent',axes1,'Color',[0 1 0],...
%     'MarkerSize',6,'LineWidth',2,'DisplayName','ref.  4')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot([loc loc],[A B],'Parent',axes1,'Color',[0 0 0],...
%     'MarkerSize',6,'LineWidth',1,'LineStyle','--',...
%     'DisplayName','time to select')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xlim(axes1,[0 D])
ylim(axes1,[A B])
% Create xlabel
xlabel('$t (day)$','FontSize',16,'FontName','Times New Roman',...
    'FontAngle','italic','Interpreter','latex');

% Create ylabel
ylabel('Pressure ($MPa$)','FontSize',16,'FontName',...
    'Times New Roman','FontAngle','italic','Interpreter','latex');

% Set the remaining axes properties
set(axes1,'FontName',...
    'Times New Roman','FontSize',14,'TickDir','both','TickLabelInterpreter',...
    'latex','XMinorTick','on','YMinorTick','on');

% Create legend
legend1 = legend(axes1,'show');
set(legend1,'Location','NorthEast','FontSize',8);
set(legend1,'Box','off');

% Print
base=[homf base_name];
%print('-djpeg90',base)
print('-depsc','-r300',base)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% NORMA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
norma=norm(ref(:,2:end))
norma=norm(ref(:,2:end)-dados(:,2:end))/norma;
fprintf('ERRO RELATIVO = %e\n',norma)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%clear
