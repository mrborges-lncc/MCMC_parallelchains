clear;
close all
loc=300;
jump=8;
N=00;
B=80.;
A=20;
% 1 m3 = 6.2898105697751 bbl
fat = 6.2898105697751;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Nch_ini = 0;
Nch_fim = 1;
Nchains = Nch_fim - Nch_ini + 1;
Nini = repmat(200, 1, Nchains);
Nfim = [400 400 600];
Nfim = Nfim(Nch_ini+1:Nch_fim+1);
Nt   = (Nfim-Nini)+1;
chains = [Nch_ini:1:Nch_fim];
base_name = 'TwoPhase3D_RW_RK';
hom  = '~/Dropbox/PROJETO_MCMC_RIGID/MCMC_parallelchains/';
homf = '~/Dropbox/PROJETO_MCMC_RIGID/paper/figuras/';
dados=load([hom 'twophaseflow/exp/pres/pres_referencia_0.dat']);
ref=dados;
home = [hom 'MonteCarlo/twophaseflow/exp/pres/'];
%
data= zeros(size(ref,1),size(ref,2),sum(Nt));
rep = [];
k = 0;
for j=1:Nchains
    n = num2str(chains(j),'%d');
    file_name =...
        [hom 'twoStage/out/nchain_' base_name n '.dat'];
    r   = load(file_name);
    rep = [rep; r(Nini(j)+1:Nfim(j)+1,2)];
    for i=Nini(j):Nfim(j)
        k = k+1;
        istr=num2str(i,5);
        fprintf('\nChain %s <=> Sample no. %s',n,istr);        
        file_name   = [home 'pres_MCMC_' base_name n '_' istr '.dat'];
        data(:,:,k) = load(file_name);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sumrep = sum(rep);
sumtot = sum(sumrep);
cont   = 0;
k      = 1;
soma   = 0.0;
soma2  = 0.0;
for j=1:Nchains
    n = 0;
    for i=Nini(j)+1:Nfim(j)+1
        dat = data(:,:,k);
        k   = k + 1;
        n   = n + 1;
        for m=1:rep(n)
            cont = cont +1;
            soma = soma + dat;
            soma2= soma2 + dat.^2;
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dmedio = soma/cont;
d2med  = soma2/cont;
erro   = sqrt(d2med - dmedio.^2);
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
    'Parent',axes1,'Color',[0.85 0.33 0.10],'MarkerSize',4,'Marker','o',...
    'LineStyle','none','DisplayName','mean','LineWidth',0.5)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plot(ref(2:end,1),ref(2:end,2),'Parent',axes1,'Color',[0.85 0.33 0.10],...
    'MarkerSize',6,'LineWidth',2,'DisplayName','ref.')
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
ylabel('Pressure ($MPa$)','FontSize',16,'FontName',...
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
base=[homf base_name '_MCMC'];
%print('-djpeg90',base)
print('-depsc','-r300',base)
pause(2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% NORMA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
norma=norm(ref(:,2:end))
norma=norm(ref(:,2:end)-dados(:,2:end))/norma;
fprintf('ERRO RELATIVO = %e\n',norma)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
