clear;
close all
loc=150;
jump=3;
N=00;
B=80.;
A=20;
% 1 m3 = 6.2898105697751 bbl
fat = 6.2898105697751;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Nch_ini = 0;
Nch_fim = 5;
Nchains = Nch_fim - Nch_ini + 1;
Nini = repmat(250, 1, Nchains);
Nfim = [530 546 570 560 553 549];
Nfim = Nfim(Nch_ini+1:Nch_fim+1)-1;
Nt   = (Nfim-Nini)+1;
chains = [Nch_ini:1:Nch_fim];
base_name = 'TwoPhase3D_RW_RK';
hom  = '~/Dropbox/PROJETO_MCMC_RIGID/MCMC_parallelchains/';
homf = '~/Dropbox/PROJETO_MCMC_RIGID/paper/figuras/';
ref  = load([hom 'MonteCarlo/twophaseflow/exp/pres/pres_ref_0.dat']);
home = [hom 'twoStage/select_prod/'];
%
data= zeros(size(ref,1),size(ref,2),sum(Nt));
rep = [];
k = 0;
for j=1:Nchains
    n = num2str(chains(j),'%d');
    file_name =...
        [hom 'twoStage/out/nchain_' base_name n '.dat']
    r   = load(file_name);
    rep = [rep; r(Nini(j)+1:Nfim(j)+1,2)];
    for i=Nini(j):Nfim(j)
        k = k+1;
        istr=num2str(i,5);
        file_name   = [home 'prod_D1_' base_name n '_' istr '.dat']
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
dados   = ref;
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
    'LineStyle','none','DisplayName','mean 1','LineWidth',0.5)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dados=ref;
%dados=load('../../MCMC/reject_prod/prod_chn0-20_16000.dat');
%dados=load('../conc/conc_amostra_0.dat');
plot(dados(2:end,1),dados(2:end,2),'Parent',axes1,'Color',[0.85 0.33 0.10],...
    'MarkerSize',6,'LineWidth',2,'DisplayName','ref.  1')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot([loc loc],[A B],'Parent',axes1,'Color',[0 0 0],...
%     'MarkerSize',6,'LineWidth',1,'LineStyle','--',...
%     'DisplayName','selection time')

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
set(legend1,'Location','NorthEast','FontSize',8);
set(legend1,'Box','off');

% Print
base=[homf 'chain_' base_name];
%print('-djpeg90',base)
print('-depsc','-r300',base)
pause(2)
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% NORMA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
norma=norm(ref(:,2:end))
norma=norm(ref(:,2:end)-dados(:,2:end))/norma;
fprintf('ERRO RELATIVO = %e\n',norma)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
