clear;
close all
loc=150;
jump=3;
N=0;
B=230;
A=70;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Nch_ini = 0;
Nch_fim = 3;
Nchains = Nch_fim - Nch_ini + 1;
Nini = repmat(1000, 1, Nchains);
Nfim = [1049 1096 1160 971];
Nfim = Nfim(Nch_ini+1:Nch_fim+1);
Nt   = (Nfim-Nini)+1;
chains = [Nch_ini:1:Nch_fim];
nome = 'TwoPhase3D_RW_RK';
nome = 'TwoPhase3D_onlyPerm_RW_RK';
base_name = ['prod_D1_' nome];
hom  = '~/Dropbox/PROJETO_MCMC_RIGID/MCMC_parallelchains/';
hom  = '~/Dropbox/PROJETO_MCMC_RIGID/MCMCrw_onlyPerm/';
homf = '~/Dropbox/PROJETO_MCMC_RIGID/paper/figuras/';
dados=load([hom 'twophaseflow/exp/pres/pres_referencia_0.dat']);
ref=dados;
home = [hom 'twoStage/select_prod/'];
%
file_name   = [home base_name '0_0.dat'];
data=load(file_name);
total=0;
for j=1:Nchains
    n = num2str(chains(j),'%d');
    pchains = load([home '../out/nchain_' nome n '.dat']);
    total = total + sum(pchains(Nini(j):Nfim(j),2));
end
data = zeros(size(data,1),size(data,2),total);
%% MEDIA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cont = 0;
for j=1:Nchains
    n = num2str(chains(j),'%d')
    pchains = load([home '../out/nchain_' nome n '.dat']);
    for i=Nini(j):Nfim(j)
        istr=num2str(i,5);
        file_name = [home base_name n '_' istr '.dat'];
        dat  = load(file_name);
        m = pchains(i,2);
        for nc = 1:m
            cont = cont + 1;
            data(:,:,cont) = dat;
        end
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dados = dmedio;
errorbar(dados(1:jump:end,1),dados(1:jump:end,2),erro(1:jump:end,2),...
    'Parent',axes1,'Color',[1 0 0],'MarkerSize',4,'Marker','o',...
    'LineStyle','none','DisplayName','mean','LineWidth',0.5)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plot(ref(2:end,1),ref(2:end,2),'Parent',axes1,'Color',[1 0 0],...
    'MarkerSize',6,'LineWidth',2,'DisplayName','ref.')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plot([loc loc],[A B],'Parent',axes1,'Color',[0 0 0],...
    'MarkerSize',6,'LineWidth',1,'LineStyle','--',...
    'DisplayName','selection time')
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
N = min(size(ref,1),size(dados,1))
norma=norm(ref(1:N,2:end))
norma=norm(ref(1:N,2:end)-dados(1:N,2:end))/norma;
fprintf('ERRO RELATIVO = %e\n',norma)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
