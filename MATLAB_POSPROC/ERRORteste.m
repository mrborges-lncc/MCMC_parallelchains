clear;
close all;
ini = 0;
fim = 79;
N   = 45000;
NV  = 1;
M   = fim-ini+1;
home= '/home/mrborges/MCMC_parRW/trunk/twoStage/';
% home= '/home/mrborges/SDUMONTruns/MCMC_parDE/trunk/twoStage/';
home= '/media/mrborges/data/MCMC_parRW/trunk/twoStage/';
home= '/media/mrborges/data/MCMC_parDE/trunk/twoStage/';
homef='/home/mrborges/MCMC_par/trunk/';
% homef='/home/mrborges/Dropbox/Eventos/2019/11th_INTERPORE/slides/';
base_name  = 'FS_RW_RK';
base_name  = 'FS_DE_RK';
% base_name  = 'blackboxRW_RK';
nome_extra = '2';
nchain=1;
my=1.250;
my0=0.0;
lwd  = 2;
xmaximo = 0;
data = [];
rep  =[];
tm   = 0;
info = [];
aux  = 0;
for i=ini:fim
    file_name =...
        [home 'error/erros_' base_name num2str(i,'%1.1d') '.dat'];
    dados=load(file_name);
    sz = size(dados,1);
    tm = [tm;sz];
    data =[data; dados];
    file_name =...
        [home 'out/nchain_' base_name num2str(i,'%1.1d') '.dat'];
    dados=load(file_name);
    dados=[dados; sz 1];
    aux = max(max(sum(dados(:,2))),aux);
    rep = [rep;dados];
    file_name =...
        [home 'in/init_stat_' base_name num2str(i,'%1.1d') '.in']
    fileID = fopen(file_name);
    [A] = fscanf(fileID,'%d %d %d %d');
    info= [info; A(1) A(3)];
    fclose(fileID);
end
tm
if(N==0)
    N = aux+aux/90;
end
% if(sum(rep(:,2))<(N-N*0.1))
%     N = sum(rep(:,2));
% end
% Create figure
figure1 = figure(1);
% fig = gcf;
% fig.Color = 'none';
% fig.InvertHardcopy = 'off';

% Create axes
axes1 = axes('Parent',figure1,'LineWidth',2,'FontSize',16,...
    'FontName','Times New Roman','FontWeight','bold',...
    'DataAspectRatio',[N/(my-my0) 3 1],'Color','none');
box(axes1,'on');
hold(axes1,'all');
mx = 0;
ny = 1e32;
final = 0;
tempo = cputime;
numcad=ini;
for i=1:M
    cor    = [(M-i)/(M-1) 0 (i-1)/(M-1)];
    cor    = [(i-1)/(M-1) (i-1)/(M-1) (i-1)/(M-1)];
    inicio = final+1;
    final  = inicio+tm(i+1)-1;
    dx     = data(inicio:final,1);
    dy     = data(inicio:final,3);
    r      = rep(inicio:final,2);
    sz     = size(r,1);
    n=0;
    if(nchain==1)
        rsum = sum(r);
        x = [1:1:rsum]';
        y = zeros(rsum,1);
    else
        x = [1:1:tm(i+1)]';
        r = 0*r+1;
    end
    for j=1:sz
        for k=1:r(j)
            n=n+1;
            x(n) = n;
            y(n) = dy(j);
        end
    end
%     mx     = max(mx,max(x));
%     my     = max(my,max(y));
%     ny     = min(ny,min(y));
    maxx = size(x,1);
    if(maxx>N)
        maxx=N
    end
    plot(x(1:maxx,1),y(1:maxx,1),'LineWidth',lwd,'Color',cor)
    numcad
    pause
    numcad=numcad+1;
    etempo = cputime-tempo
end
% Create xlabel
xlabel('Accepted iterations','FontWeight','bold','Interpreter','latex',...
    'LineWidth',lwd,'FontSize',16,...
    'FontName','Times New Roman');
% Create ylabel
ylabel('$\mathsf{Er}$','Interpreter','latex','LineWidth',lwd,...
    'FontSize',16,'FontName','Times New Roman');
%ylim([ny*0.9 my]);
ylim([my0 my]);
xlim([0 maxx]);
set(gcf,'PaperPositionMode','auto');
base=[homef 'figuras/error1_' base_name '_' nome_extra]
print('-depsc','-r300',base);
%print('-djpeg','-r300',base);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(NV>1)
    my=1.250;
    % Create figure
    figure2= figure;
    fig = gcf;
    fig.Color = 'none';
    fig.InvertHardcopy = 'off';
% Create axes
    axes2 = axes('Parent',figure2,'LineWidth',2,'FontSize',16,...
        'FontName','Times New Roman','FontWeight','bold',...
        'DataAspectRatio',[N/(my-my0) 3 1],'Color','none');
    box(axes2,'on');
    hold(axes2,'all');

    mx = 0;
    ny =1e32;
    final = 0;
    for i=1:M
        cor    = [(M-i)/(M-1) 0 (i-1)/(M-1)];
        cor    = [(i-1)/(M-1) (i-1)/(M-1) (i-1)/(M-1)];
        inicio = final+1;
        final  = inicio+tm(i+1)-1;
        dx      = data(inicio:final,1);
        dy      = data(inicio:final,4);
        r      = rep(inicio:final,2);
        sz     = size(r,1);
        x=[];
        y=[];
        n=0;
        for j=1:sz
            if(nchain~=1)r(j)=1;end
            for k=1:r(j)
                n=n+1;
                x = [x;n];
                y = [y;dy(j)];
            end
        end
    %     mx     = max(mx,max(x));
    %     my     = max(my,max(y));
    %     ny     = min(ny,min(y));
        plot(x,y,'LineWidth',lwd,'Color',cor)
        hold on
    end
    % Create xlabel
    xlabel('Accepted iterations','FontWeight','bold',...
        'Interpreter','latex','LineWidth',lwd,'FontSize',16,...
        'FontName','Times New Roman');
    % Create ylabel
    ylabel('$\mathsf{E}_{\ 2}$','Interpreter','latex',...
        'LineWidth',lwd,'FontSize',16,...
        'FontName','Times New Roman');
    %ylim([ny*0.9 my*1.1]);
    ylim([my0 my]);
    xlim([0 N]);
    base=[homef 'figuras/error2_' base_name '_' nome_extra]
    set(gcf,'PaperPositionMode','auto');
    print('-depsc','-r300',base);
%     print('-djpeg','-r300',base);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(NV>2)
    % Create figure
    figure3= figure;
    fig = gcf;
    fig.Color = 'none';
    fig.InvertHardcopy = 'off';    % Create axes
    axes3 = axes('Parent',figure3,'LineWidth',2,'FontSize',16,...
        'FontName','Times New Roman','FontWeight','bold',...
        'DataAspectRatio',[N/(my-my0) 3 1],'Color','none');
    box(axes3,'on');
    hold(axes3,'all');

    mx = 0;
    ny = 1e32;
    final = 0;
    for i=1:M
        cor    = [(M-i)/(M-1) 0 (i-1)/(M-1)];
        cor    = [(i-1)/(M-1) (i-1)/(M-1) (i-1)/(M-1)];
        inicio = final+1;
        final  = inicio+tm(i+1)-1;
        dx     = data(inicio:final,1);
        dy     = data(inicio:final,5);
        r      = rep(inicio:final,2);
        sz     = size(r,1);
        x=[];
        y=[];
        n=0;
        for j=1:sz
            if(nchain~=1)r(j)=1;end
            for k=1:r(j)
                n=n+1;
                x = [x;n];
                y = [y;dy(j)];
            end
        end
    %     mx     = max(mx,max(x));
    %     my     = max(my,max(y));
    %     ny     = min(ny,min(y));
        plot(x,y,'LineWidth',lwd,'Color',cor)
        hold on
    end
    % Create xlabel
    xlabel('Accepted iterations','FontWeight','bold',...
        'Interpreter','latex',...
        'LineWidth',lwd,'FontSize',16,...
        'FontName','Times New Roman');
    % Create ylabel
    ylabel('$\mathsf{E}_{\ 3}$','Interpreter','latex',...
        'LineWidth',lwd,'FontSize',16,...
        'FontName','Times New Roman');
    %ylim([ny*0.9 my*1.1]);
    ylim([my0 my]);
    xlim([0 N]);
    base=[homef 'figuras/error3_' base_name '_' nome_extra]
    %pause
    set(gcf,'PaperPositionMode','auto');
    print('-depsc','-r300',base);
%     print('-djpeg','-r300',base);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create figure
figure4= figure;
    fig = gcf;
    fig.Color = 'none';
    fig.InvertHardcopy = 'off';
N=max(tm)+1;
my=1.25;
% Create axes
axes3 = axes('Parent',figure4,'LineWidth',2,'FontSize',16,...
    'FontName','Times New Roman','FontWeight','bold',...
    'DataAspectRatio',[N/(my-my0) 3 1],'Color','none');
box(axes3,'on');
hold(axes3,'all');

mx = 0;
ny = 1e32;
final = 0;
nchain = 10;
Ninicial = 100;
Y=[];
%N=max(tm)
for i=1:M
    cor    = [(M-i)/(M-1) 0 (i-1)/(M-1)];
    cor    = [(i-1)/(M-1) (i-1)/(M-1) (i-1)/(M-1)];
    inicio = final+1;
    final  = inicio+tm(i+1)-1;
    dx     = data(inicio:final,1);
    dy     = data(inicio:final,2);
    r      = rep(inicio:final,2);
    sz     = size(r,1);
    x=[];
    y=[];
    n=0;
    for j=1:sz
        if(nchain~=1)r(j)=1;end
        for k=1:r(j)
            n=n+1;
            x = [x;n];
            y = [y;dy(j)];
        end
    end
%     mx     = max(mx,max(x));
%     my     = max(my,max(y));
%     ny     = min(ny,min(y));
    plot(x,y,'LineWidth',lwd,'Color',cor)
    hold on
    Y = [Y; y(Ninicial:end,1)];
end
% Create xlabel
xlabel('Accepted iterations','FontWeight','bold',...
    'Interpreter','latex',...
    'LineWidth',lwd,'FontSize',16,...
    'FontName','Times New Roman');
% Create ylabel
ylabel('$\mathsf{Er}$','Interpreter','latex',...
    'LineWidth',lwd,'FontSize',16,...
    'FontName','Times New Roman');
%ylim([ny*0.9 my]);
ylim([my0 my]);
xlim([0 N]);
base=[homef 'figuras/error4_' base_name '_' nome_extra]
%pause
set(gcf,'PaperPositionMode','auto');
print('-depsc','-r300',base);
% print('-djpeg','-r300',base);
tm
taxa=100.0*double(info(:,2))./double(info(:,1));
fprintf('#################################\n');
fprintf('#################################\n');
fprintf('TAXA DE ACEITACCAO MÉDIA: %1.2f\n',mean(taxa))
fprintf('#################################\n');
fprintf('#################################\n');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
est=10;
errm = zeros(M,1);
nchain=1;
if(est==1)
    Ninicial = 1;
    Y=[];
    final=0;
    errm = zeros(M,1);
    %N=max(tm)
    for i=1:M
        cor    = [(M-i)/(M-1) 0 (i-1)/(M-1)];
        cor    = [(i-1)/(M-1) (i-1)/(M-1) (i-1)/(M-1)];
        inicio = final+1;
        final  = inicio+tm(i+1)-1;
        dx     = data(inicio:final,1);
        dy     = data(inicio:final,2);
        r      = rep(inicio:final,2);
        sz     = size(r,1);
        x=[];
        y=[];
        n=0;
        for j=1:sz
            if(nchain~=1)r(j)=1;end
            for k=1:r(j)
                n=n+1;
                x = [x;n];
                y = [y;dy(j)];
            end
        end
        Ninicial = floor(size(y,1)/2);
        errm(i,1) = mean(log(y(Ninicial:end,1)));
        %Y = [Y; y(Ninicial:end,1)];
    end
end
% Create figure
figure1 = figure(10);
    fig = gcf;
    fig.Color = 'none';
    fig.InvertHardcopy = 'off';

% Create axes
axes3 = axes('Parent',figure1,'LineWidth',2,'XTickLabel','',...
    'XTick',zeros(1,0),'FontSize',12,'FontName','Times New Roman',...
    'Position',[0.13 0.1529 0.6 0.8]);
boxplot(axes3,errm,'BoxStyle','outline','symbol','+',...
    'colors',[0 0 0],'widths',0.85)
ylabel('mean log error','FontWeight','bold',...
    'FontSize',14,...
    'FontName','Times New Roman');
base=[homef 'figuras/Boxplot_' base_name '_' nome_extra];
set(gcf,'PaperPositionMode','auto');
print('-depsc','-r300',base);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
q3 = quantile(errm,0.750);
median(errm);
tol = q3+2.0*iqr(errm);
% tol = median(errm)+2.0*iqr(errm);
% hold on
% y=[tol tol];
% x=[0.0 2];
% plot(axes3,x,y,'LineWidth',lwd,'LineStyle','-','Color',[1 0 0])
% ylabel('mean log error','FontWeight','bold',...
%     'FontSize',14,...
%     'FontName','Times New Roman');
outliers=[];
for i=1:M
    if(errm(i,1)>tol)
        outliers = [outliers; i-1];
    end
end
outliers;
MIN = min(tm(2:end))
clear
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%clear
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%clear
