clear;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PARAMETROS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
speedshock=3.6475;
sstar     =0.3854;
poro      =0.104495;
%
% speedshock=1;
% sstar     =0.3854;
% poro      =0.1077;
%
sl=0.21;
auxd=1.0;
dt=10;          % intervalo de tempo
Lx=1.0;
nx=200;
dx=Lx/nx;
h=dx*0.5;
TEMPO_CHEGADA = poro*Lx/speedshock;
PROD=100.7;
PROD_MIN=0.0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FIGURAS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
name_fig = '../figuras/prod_het';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('----------------')
disp(' LOADING OUTPUT DATA ');
%%%%%% NUMERO DE ARQUIVOS PARA COMPARACAO e LEITURA %%%%%%%%%%%%%%%%%%%%%%%
N=3;
pr=[];
nmin=1e+10;
nmax=-1e10;
name=[];
nl=0;
for i=1:N
    [FILENAME,PATHNAME] = uigetfile('../prod/prod_h*.dat', 'LOAD DATA');
    file=[sprintf('%s%s',PATHNAME,FILENAME)];
    name=[name; file]
    p = load(file);
    nmin = min(nmin,length(p));
    nmax = max(nmax,length(p));
    pr = [pr; p];
    nl= [nl; length(pr)]
    clear p
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
prod=zeros(N,nmax,2);
for i=1:N
    prod(i,:,:)=pr(nl(i)+1:nl(i+1),:);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Create figure
figure1 = figure('XVisual',...
    '0x3b (TrueColor, depth 32, RGB mask 0xff0000 0xff00 0x00ff)');

% Create axes
y_ini=min(min(prod(:,:,2)));
y_lim=max(max(prod(:,:,2)))*1.1;
x_ini=min(min(prod(:,:,1)));
x_lim=max(max(prod(:,:,1)));

axes1 = axes(...
  'DataAspectRatio',[1/(y_lim-y_ini) 2/(x_lim-x_ini) 1],...
  'Parent',figure1,'FontSize',16,'FontName','Times New Roman');
ylim(axes1,[y_ini y_lim]);
xlim(axes1,[x_ini x_lim]);
box('on');
grid('on');
hold('all');

for i=1:N
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if(i==1)
        simbol='o'
        if(name(i,end-6)=='C')
            leg='$KC$'
        end
        if(name(i,end-6)=='P')
            leg='$\phi$'
        end
        if(name(i,end-6)=='K')
            leg='$k$'
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if(i==2)
        simbol='+'
        if(name(i,end-6)=='C')
            leg='$KC$'
        end
        if(name(i,end-6)=='P')
            leg='$\phi$'
        end
        if(name(i,end-6)=='K')
            leg='$k$'
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if(i>2)
        simbol='s'
        if(name(i,end-6)=='C')
            leg='$KC$'
        end
        if(name(i,end-6)=='P')
            leg='$\phi$'
        end
        if(name(i,end-6)=='K')
            leg='$k$'
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    plot(prod(i,:,1),prod(i,:,2),'Parent',axes1,...
        'Marker',simbol,'DisplayName',leg);
end
% Create xlabel
xlabel('$t$','Interpreter','latex','FontSize',24,...
    'FontName','Times New Roman');

% Create ylabel
ylabel('Prod. oleo','Interpreter','latex','FontSize',22,...
    'FontName','Times New Roman');

% Create legend
legend1 = legend(axes1,'show');
set(legend1,'Position',[0.75 0.625 0.07539 0.07296],...
    'Interpreter','latex');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
name_fig = [name_fig];
print('-depsc',name_fig);
