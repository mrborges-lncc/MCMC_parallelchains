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
PROD=0.5;
PROD_MIN=0.0;
jump=5;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FIGURAS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
name_fig = '../figuras/prod_het';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('----------------')
disp(' LOADING OUTPUT DATA ');
%%%%%% NUMERO DE ARQUIVOS PARA COMPARACAO e LEITURA %%%%%%%%%%%%%%%%%%%%%%%
N=5;
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
    nl= [nl; length(pr)];
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
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     if(name(i,end-6)=='C')
%         leg='KC';
%     end
%     if(name(i,end-6)=='P')
%         leg='\phi';
%     end
%     if(name(i,end-6)=='K')
%         leg='k';
%     end
%     if(name(i,end-6)=='d')
%         leg='k, \phi uncor.';
%     end
%     if(name(i,end-6)=='g')
%         leg='hom.';
%     end
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if(name(i,end-6)=='C')
        leg='case iv';
    end
    if(name(i,end-6)=='P')
        leg='case i';
    end
    if(name(i,end-6)=='K')
        leg='case ii';
    end
    if(name(i,end-6)=='d')
        leg='case iii';
    end
    if(name(i,end-6)=='g')
        leg='hom.';
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if(i==1)
        simbol='o';
        cor=[0 0 1];
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if(i==2)
        simbol='+';
        cor=[1 0 0];
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if(i==3)
        simbol='s';
        cor=[0 1 0];
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if(i==4)
        simbol='^';
        cor=[0.07843 0.1686 0.549];
        cor=[0 0 0];
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if(i>4)
        simbol='none';
        cor=[0 0 0];
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    plot(prod(i,1:jump:end,1),prod(i,1:jump:end,2),'Parent',axes1,...
        'Color',cor,'Marker',simbol,'MarkerSize',5,'DisplayName',leg);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
% Create xlabel
xlabel('$t$','Interpreter','latex','FontSize',24,...
    'FontName','Times New Roman');

% Create ylabel
ylabel('Oil prod.','Interpreter','tex','FontSize',22,...
    'FontName','Times New Roman');

% Create legend
legend1 = legend(axes1,'show');
set(legend1,'Position',[0.775 0.60 0.07539 0.07296],...
    'Interpreter','tex','FontAngle','italic');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
name_fig = [name_fig];
print('-depsc',name_fig);
