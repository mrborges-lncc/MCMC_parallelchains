clear;
close all;
ini = 0;
fim = 3;
N   = 0;
NV  = 2;
M   = fim-ini+1;
home= '/home/mrborges/MCMCrw/twoStage/';
% home= '~/Dropbox/PROJETO_MCMC_RIGID/MCMC_parallelchains/twoStage/';
% home= '../twoStage/';
%
homef='~/Dropbox/PROJETO_MCMC_RIGID/paper/figuras/';
% homef='~/MCMC_parallelchains/';
% homef='../';
base_name = 'TwoPhase3D_RW_RK';
nome_extra = '';
nchain=1;
razao = 4;
my  = [0.02; 0.08; 0.3];
my0 = [0.00; 0.00; 0.00];
lwd = 2;
xmaximo= 0;
transp = 0;
data = [];
rep  = [];
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

if(N==0)
    N = ceil(aux+aux/90);
end

nvar = size(data,2) - 2;

% Create figure
for nf = 1:nvar + 1
    figure1 = figure('PaperOrientation','landscape','PaperSize',[11 5]);
    if(transp==1)
        fig = gcf;
        fig.Color = 'none';
        fig.InvertHardcopy = 'off';
    end
    % Create axes
    dasp = [N/(my(nf)-my0(nf)) razao 1];
    axes1 = axes('Parent',figure1,'LineWidth',1,'FontSize',16,...
        'FontName','Times New Roman','FontWeight','bold',...
        'DataAspectRatio',dasp,'Color','none');
    % Set the remaining axes properties
    set(axes1,'LineWidth',1,'TickDir','both',...
        'TickLabelInterpreter','latex','XMinorTick','on',...
        'YMinorTick','on');
    box(axes1,'on');
    hold(axes1,'all');
    final = 0;
    for i=1:M
        cor    = [(M-i)/(M-1) 0 (i-1)/(M-1)];
        %cor    = [(i-1)/(M-1) (i-1)/(M-1) (i-1)/(M-1)];
        inicio = final+1;
        final  = inicio+tm(i+1)-1;
        dx     = data(inicio:final,1);
        dy     = data(inicio:final,nf+1);
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
        maxx = size(x,1);
        if(maxx>N)
            maxx=N;
        end
        plot(x(1:maxx,1),y(1:maxx,1),'LineWidth',lwd,'Color',cor)
    end
    % Create xlabel
    xlabel('Accepted iterations','FontWeight','bold','Interpreter','latex',...
        'LineWidth',lwd,'FontSize',16,...
        'FontName','Times New Roman');
    % Create ylabel
    name = ['$\mathsf{E}_{ ' num2str(nf-1,'%d') '}$'];
    if nf == 1, name = ['$\mathsf{Er}$']; end
    ylabel(name,'Interpreter','latex','LineWidth',lwd,...
        'FontSize',16,'FontName','Times New Roman');

    ylim([my0(nf) my(nf)]);
    xlim([0 N]);
    set(gcf,'PaperPositionMode','manual','PaperPosition',[0.25 0.25 2.25*razao 3]);
    base=[homef 'error' num2str(nf-1,'%d') '_' base_name nome_extra]
    print('-depsc','-r300',base);
    pause(1); clf; close all;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%clear
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%clear
tm
taxa=100.0*double(info(:,2))./double(info(:,1));
fprintf('#################################\n');
fprintf('#################################\n');
fprintf('TAXA DE ACEITACCAO MÉDIA: %1.2f\n',mean(taxa))
fprintf('#################################\n');
fprintf('#################################\n');
clear
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%clear
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%clear
