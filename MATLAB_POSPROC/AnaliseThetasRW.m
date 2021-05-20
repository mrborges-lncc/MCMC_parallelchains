close all;
clear;
analysis = 1;
convanalysis = 1;
tempo = cputime;
d  = 40; % dimension of multivariate Gausesian
Nc0= 0; % Numero da cadeia inicial
Nc = 79; % Numero da cadeia final
NT = Nc-Nc0+1; % Numero total de cadeias
NCAD = [Nc0:1:Nc];
MM   = 200;
home = '/home/mrborges/MCMC_parDE/trunk/';
home= '/media/mrborges/data/MCMC_parDE/trunk/';
home= '/media/mrborges/data/MCMC_parRW/trunk/';
home= '/home/mrborges/MCMC_parRW/trunk/';
homef= '/home/mrborges/MCMC_par/trunk/';
homef='/home/mrborges/Dropbox/Eventos/2019/11th_INTERPORE/slides/';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
data = [];
tm = [];
rep=[];
base_name = 'FS_DE_RK';
base_name = 'FS_RW_RK';
nome_extra = '';
refile = [home 'simuladorRigido/exp/conc/conc_ref_'];
for i=Nc0:Nc
    file_name = ...
        [home 'twoStage/error/erros_' base_name num2str(i,'%1.1d') '.dat'];
    dados=load(file_name);
    sz = size(dados,1);
    tm = [tm;sz];
    Ncurv = size(dados,2)-2;
    file_name =...
        [home 'twoStage/out/nchain_' base_name num2str(i,'%1.1d') '.dat'];
    dados=load(file_name);
    dados=[dados; sz 1];
    rep = [rep;dados];
end
tm=tm-1;
Nfim = tm-1;
% Nfim = Nfim*0+202;
Nini =0*Nfim+MM;
N=max(tm)+1;
clear dados
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% DADOS DE REFERENCIA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
szref  = 1;
ref    = [];
for i=1:Ncurv
    nomef = [refile num2str(i-1,3) '.dat'];
    dataref = load(nomef);
    sz = size(dataref,1);
    szref = max(szref,sz);
    ref = [ref dataref];
end
szref = size(ref)-1;
ref = ref(1:szref,:);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
variavel  = [home 'twoStage/select_thetas/theta_v'];
variavel1 = [home 'twoStage/out/nchain_'];
variavel2 = [home 'twoStage/select_prod/prod_D'];
nchain=zeros(N,NT);
szz = size(ref,2);
sz  = size(ref,1);
nk=0;
SOMA = [];
mu  = zeros(1,d);
mu2 = zeros(1,d);
scur= zeros(sz,szz);
dcur= zeros(sz,szz);
nk = 0;
for i=1:NT
    nfile =...
        [variavel1 base_name num2str(NCAD(i),3) '.dat'];
    rep = load(nfile);
    nchain = rep(:,2);
    for k=Nini(i):Nfim(i)
        nfile =...
            [variavel '1_' base_name num2str(NCAD(i),3) '_' num2str(k,5) '.dat']
        nk=nk+nchain(k);
        x = load(nfile)';
        nfile =...
            [variavel2 '1_' base_name num2str(NCAD(i),3) '_' num2str(k,5) '.dat'];
        curva = load(nfile);
        curva = curva(1:szref(1),:);
        scur  = scur+nchain(k)*curva;
        dcur  = dcur+nchain(k)*curva.^2;
        mu    = mu + nchain(k)*x;
        mu2   = mu2 + nchain(k)*x.^2;
    end
end
mu  = mu/nk;
mu2 = mu2/nk;
scur= scur/nk;
dcur= dcur/nk;
sd  = sqrt(dcur(:,2:end)-scur(:,2:end).^2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
info=[];
for i=Nc0:Nc
    file_name = [home 'twoStage/error/erros_' base_name num2str(i,'%1.1d') '.dat'];
    dados=load(file_name);
    sz = size(dados,1);
    tm = [tm;sz];
%     data =[data; dados];
    file_name = [home 'twoStage/out/nchain_' base_name num2str(i,'%1.1d') '.dat'];
    dados=load(file_name);
    dados=[dados; sz 1];
    rep = [rep;dados];
    file_name = [home 'twoStage/in/init_stat_' base_name num2str(i,'%1.1d') '.in']
    fileID = fopen(file_name);
    [A] = fscanf(fileID,'%d %d %d %d');
    info= [info; A(1) A(3)];
    fclose(fileID);
end
taxa=100.0*double(info(:,2))./double(info(:,1));
AR = mean(taxa);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xdata = scur(:,1);
B = max(max(ref(:,2:end)))*1.1;
A = min(min(ref(:,2:end)));
A = 0.;
B = 1.25;
C = min(xdata)*1.;
D = max(xdata)*1.;
C = 0.0;
D = 300;
% Create figure
figure1 = figure(1);

% Create axes
axes1 = axes('Parent',figure1,'LineWidth',2,'FontSize',18,...
    'FontName','Times New Roman','FontWeight','bold',...
    'DataAspectRatio',[1 3*(B-A)/(D-C) 1],...
    'Color','none');
box(axes1,'on');
hold(axes1,'all');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pula=1;
for j=1:NT
    for k=Nini(j):Nfim(j)
        nfile =...
            [variavel2 '1_' base_name num2str(NCAD(j),3) '_' num2str(k,5) '.dat']
        curva = load(nfile);
%         for i=2:szz
%             plot(curva(:,1),curva(:,i),'-','LineWidth',2,'Color',[0.85 0.85 0.85]);
%         end
    end
end
jump = 10;
for i=2:pula:szz
    errorbar(scur(1:jump:end,1),scur(1:jump:end,i),sd(1:jump:end,i-1),...
        'Marker','o','MarkerSize',8,...
        'LineWidth',2,'Color',[1-i/szz 0 i/szz],'LineStyle','none')
end
for i=2:pula:szz
    plot(ref(:,1),ref(:,i),'-','LineWidth',3,'Color',[1-i/szz 0 i/szz]);
end
% Uncomment the following line to preserve the X-limits of the axes
xlim(axes1,[min(xdata)*1. max(xdata)*1.1]);
% Uncomment the following line to preserve the Y-limits of the axes
ylim(axes1,[A B]);
xlim(axes1,[C D]);
box(axes1,'on');
hold(axes1,'all');
% Create xlabel
xlabel('x','FontSize',22,'FontName','Times New Roman','FontAngle','italic');
% Create ylabel
ylabel('$s_{\! w}$','FontSize',24,'FontName','Times New Roman',...
    'Interpreter','latex','FontAngle','italic');
% Create textbox
local = [0.22 0.15 0.18 0.067];
local = [0.65 0.8 0.18 0.067];
local = [0.2 0.65 0.18 0.067];
local = [0.165 0.625 0.25 0.067];
annotation(figure1,'textbox',...
    local,...
    'String',{['AR = ' num2str(AR,'%4.1f') ' %']},...
    'HorizontalAlignment','center',...
    'FontSize',14,...
    'FontName','Times',...
    'FitBoxToText','off',...
    'LineStyle','none');
nome = [homef 'figuras/MCMCpost_'];
name = [nome base_name '_' nome_extra]
fig = gcf;
fig.Color = 'none';
fig.InvertHardcopy = 'off';
set(gcf,'PaperPositionMode','auto');
print('-depsc','-r300',name);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(convanalysis == 1)
    N = max(info(:,1));
    x = zeros(1,d);
    for i=Nc0:Nc
        file_name = ...
            [home 'twoStage/error/erros_' base_name num2str(i,'%1.1d') '.dat'];
        dados=load(file_name);
        sz = size(dados,1);
        tm = [tm;sz];
        Ncurv = size(dados,2)-2;
        file_name =...
            [home 'twoStage/out/nchain_' base_name num2str(i,'%1.1d') '.dat'];
        dados=load(file_name);
        dados=[dados; sz 1];
        rep  = [rep;dados];
    end
    for nd=1:Ncurv
        X = zeros(N,d,NT);
        for i=1:NT
            nel=0;
            nfile =...
                [variavel1 base_name num2str(NCAD(i),3) '.dat']
            rep = load(nfile);
            nchain = rep(:,2);
            nchain = [nchain; N-sum(nchain)];
            for k=0*Nini(i):Nfim(i)+1
                nfile =...
                    [variavel num2str(nd,3) '_' base_name...
                    num2str(NCAD(i),3) '_' num2str(k,5) '.dat'];
                x = load(nfile)';
                for nrep=1:nchain(k+1)
                    nel = nel+1;
                    X(nel,:,i) = x;
                end
            end
        end
%% convergencia %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        cv1 = [];
        xcv = [];
        V1  = [];
        W1  = [];
        b = N/20;
        for i=b:b:N
            ini = i/2;
            fim = i;
            xcv = [xcv; fim];
%             [cv,v,w] = convergeR(X(ini:fim,:,:));
            [cv,v,w] = convergeRMatrix(X(ini:fim,:,:));
            cv1 = [cv1; cv];
            W1  = [W1; w];
            V1  = [V1; v];
        end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Rfigure(xcv,cv1)
        nome = [homef 'figuras/MCMC_R_'];
        name = [nome base_name '_' nome_extra 'v' num2str(nd,3)]
        set(gcf,'PaperPositionMode','auto');
        print('-depsc','-r300',name)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Rfigure2(xcv,(W1),(V1))
        nome = [homef 'figuras/MCMC_WV_'];
        name = [nome base_name '_' nome_extra 'v' num2str(nd,3)]
        set(gcf,'PaperPositionMode','auto');
        print('-depsc','-r300',name)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
end
% fprintf('\n Desvios médios: sigma_1 = %f; sigma_2 = %f; sigma_3 = %f \n',...
%     mean(sd1),mean(sd2),mean(sd3))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%clear