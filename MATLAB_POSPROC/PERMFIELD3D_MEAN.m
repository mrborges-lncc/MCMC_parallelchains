%permfield.m
% -----------
escala=1;
clear;
%prompt2={'Diretorio atual figuras/: '};
%Answers2=inputdlg(prompt2,'Nome base para os arquivos',1);
%base_name = char(Answers2);
%base_aux = '../figuras/';
%base_aux = '/home/mrborges/Congressos/cilamce2011/figuras/';
%base_aux = '/home/mrborges/THERMOHIDRO/simulador/figuras/';
base_aux = '../figuras/'
N=100;
NFi=99;
NFf=149;
Lz=2000;
nz=100;
dz=Lz/N;
imprK=10; % se igual a 1 gera o grafico de K
stat=10;  % se igual a 1 gera os histogramas das vas.
tipo=1 % 1 se perm e 2 se phi
permfigs=1; % se igual a 1 gera os graficos 2D
TIPOINPUT = 1; % ==1 entrada por arquivo
TIPOINPUT2 = 10; % ==1 entrada por arquivo
filme=10;
mi=-3;
ma=3;
controle=0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
M=4.7887e-14;
S=1.0;
%M=1.0;
% S=0.2918;
% M=0.2396;
%S=1.0
%M=1.0
% M=-1.3671
% S=0.0791
%S=1./sqrt(1.2066)*1.0
%
if(tipo==2)
    grafnorm=0;
    strtipo='\phi';
end
if(tipo==1)
    grafnorm=1;
    strtipo='{k}';
    %strtipo='${E}$';
end
%**************************************************************************
%** ENTRADA DE DADOS DO CONDICIONAMENTO ***********************************
% vet(n,i) coordenada i da posicao do dado condicionado n
inp = load('../gera_KL/MATLAB/input_cond.dat');
%
if(TIPOINPUT2 == 1)
    name = '../simul_comp/exp/input_cond.in'
%     name = '../simuladorBTMM/exp/input_cond.in'; 
    fd = fopen(name,'r');
    np=fscanf(fd,'%d');
    data=zeros(np,3);
    aux=fscanf(fd,'%f');
    k=0;
    for i=1:np
        for j=1:3
            k=k+1;
            data(i,j)=aux(k);
        end
    end
    vet2 = data(:,1:2)
end
if(TIPOINPUT == 1)
    vet = inp(:,1:2)
    dados=inp(:,3);
    clear inp
else
    vet=[];
end
vet=vet+1e-6;
dados=zeros(size(vet,1),1);
pnode=zeros(size(vet,1),1);

disp('----------------')
disp(' LOADING FIELD ');

% [FILENAME, PATHNAME] =uigetfile('../forecast/3Dfields/fields/*.dat', 'LOAD DATA');
% [FILENAME, PATHNAME] =uigetfile('../simul_comp/exp/fields/*.dat', 'LOAD DATA');
% [FILENAME, PATHNAME] =uigetfile('../simuladorBTMM/exp/fields/*.dat', 'LOAD DATA');
% [FILENAME, PATHNAME] =uigetfile('../twoStage/select_fields/*.dat', 'LOAD DATA');
% [FILENAME, PATHNAME] =uigetfile('../gera_KL/MATLAB/campos/*.dat', 'LOAD DATA');
% line_file=sprintf('%s%s', PATHNAME,FILENAME);
line_file='../forecast/3Dfields/fields/MC_0_0.dat'
line_file='../forecast/FORTRAN_LBD3D/campos/MCMC_0_0.dat'
n=0;
%
% for i=size(line_file,2):-1:1
%    a=line_file(i);
%    if(a=='.')
%        in=line_file(i-1);
%        for j=i-2:-1:1
%            a=line_file(j);
%            if(a=='_')
%                str_k=j;
%                i=0;
%                break
%            else
%                in=[a in]
%            end
%        end
%    end
% end
sgn=0;
for i=size(line_file,2):-1:1
    a=line_file(i);
    if(a=='_')
        sgn=sgn+1;
        line_file(1:i)
    end
    if(sgn==2)
        str_k=i
        break
    end
end

%str_k=34
file_base=line_file(1:str_k)
%in='0'
%
for i=size(file_base,2):-1:1
    a=file_base(i);
    if(a=='/')
        str_k=i+1;
        break
    end
end
base_name = file_base(str_k:end)
%
ini=0;
fim=ini+N-1;
for nf=NFi:NFf
    cont=0;
    for II=ini:1:fim
        line_file = [file_base num2str(nf,5) '_' num2str(II,5) '.dat']
        fid = fopen(line_file,'r');
        mattamp = fscanf(fid,'%f');
        cont=cont+1;
        disp('file loaded.')
        fclose(fid);
        inf = mattamp(1:4);
        Lx = inf(1);
        Ly = inf(2);
        nx = inf(3);
        ny = inf(4);
        dx = Lx/nx;
        dy = Ly/ny;
        if(nf==NFi&&II==ini)
            permmap =zeros(ny,nx,N);
            permmapv=zeros(ny,nx,N);
        end
        if(abs(nx-ny)<1e-8)
            quad=1;
        else
            if(nx/ny==2)
                quad=0;
            else
                quad=2;
            end
        end
        mattamp = mattamp(9:length(mattamp));
        k=0;
        for j=ny:-1:1
            k=k+1;
            if(mattamp(k)~=ny-j)
                disp('erro1')
                break
            end
            for i=1:nx
                k=k+1;
                permmap(j,i,cont) =permmap(j,i,cont)+(S*mattamp(k));
                permmapv(j,i,cont)=permmapv(j,i,cont)+(S*mattamp(k))^2;
    %             permmap(j,i,cont)=M*exp(S*mattamp(k));
            end
            k=k+1;
            if(mattamp(k)~=192837465)
                disp('erro2')
                break
            end       
        end
        clear mattamp
        s = size(permmap);
        x = s(:,2);
        y = s(:,1);
        if(controle==0)
            ma=max(max(permmap))*1.01;
            mi=min(min(permmap))*1.01;
        end
%         if(quad==1)
%             posi=[0.02 0.1 0.8073 0.8026];
%             posi2=[0.75 0.1 0.045 0.8];
%         end
%         if(quad==0)
%             posi=[0.0976 0.3 0.76 0.8026];
%      %       posi=[0.09766 0.3 0.8073 0.8026];
%             posi2=[0.87 0.49 0.045 0.425];
%         end
%         if(quad==2)
%             posi=[0.09766 0.3 0.8073 0.8026];
%             posi2=[0.925 0.565 0.045 0.27];
%         end
% 
%         media=mean(mean(mean(permmap)))
%         vd=reshape(permmap,nx*ny*N,1);
%         variancia=var(vd)
%         std=sqrt(variancia)
        clear vd mattamp
    end
end
permmap =permmap/(NFf-NFi+1);
permmapv=permmapv/(NFf-NFi+1)-permmap.^2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x=dx:dx:Lx;
y=Ly-dy:-dy:0;
% x=Lx-dx:-dx:0;
% y=dy:dy:Ly;
z=dz:dz:Lz;
[X,Y,Z]=meshgrid(x,y,z);
xslice = [dx, Lx/2, Lx-dx];
yslice = [dy];
zslice = [25*dz, Lz/2, Lz-dz];
xslice = [dx,Lx/2 Lx-dx];
yslice = [Ly/2];
zslice = [dz Lz/4 Lz/2, 3*Lz/4, Lz-dz];
% Create figure
figure1 = figure('XVisual',...
    '0xcc (TrueColor, depth 24, RGB mask 0xff0000 0xff00 0x00ff)');

% Create axes
% Create axes
axes1 = axes('Parent',figure1,'XDir','reverse','FontSize',8,...
    'FontName','Times New Roman',...
    'DataAspectRatio',[1 1 3],...
    'CLim',[-2.5 2.5],...
    'CameraViewAngle',12.5);
%view(axes1,[-168.5 -69]);
%view(axes1,[-164.5 -37]);
%view(axes1,[150 -45]);
v=[150 -60];
view(axes1,v);
box(axes1,'on');
grid(axes1,'on');
hold(axes1,'all');

slice(X,Y,Z,permmap,xslice,yslice,zslice), shading flat

% Create xlabel
xlabel('$x$','Interpreter','latex','HorizontalAlignment','right',...
    'FontSize',10,...
    'FontName','Times New Roman');

% Create ylabel
ylabel('$y$','Interpreter','latex','FontSize',10,...
    'FontName','Times New Roman',...
    'HorizontalAlignment','left');

% Create zlabel
zlabel('Time',...%'Interpreter','latex',...%
    'FontSize',10,...
    'FontName','Times New Roman',...
    'FontAngle','italic');
% Create colorbar
colorbar('peer',axes1,...
    [0.58 0.35 0.01 0.35],'FontSize',8,...
    'FontName','times');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    posi2=[0.75 0.27 0.045 0.5];
    if(TIPOINPUT==1)
        plot3(vet(:,1),vet(:,2),0.0*vet(:,1),...
            'Parent',axes1,'MarkerFaceColor',[0 0 0],...
            'MarkerEdgeColor',[0 0 0],...
            'MarkerSize',6,...
            'Marker','+','LineStyle','none');
% %         %% Create colorbar
%         colorbar1 = colorbar('peer',...
%           axes1,posi2,...
%           'Box','on',...
%           'FontName','times',...
%           'FontSize',12);
    end
    if(TIPOINPUT2==1)
        plot3(vet2(:,1),vet2(:,2),0.0*vet2(:,1),...
            'Parent',axes1,'MarkerFaceColor',[0 0 0],...
            'MarkerEdgeColor',[1 1 1],...
            'MarkerSize',4,...
            'Marker','o','LineStyle','none');
%         %% Create colorbar
%         colorbar1 = colorbar('peer',...
%           axes1,posi2,...
%           'Box','on',...
%           'FontName','times',...
%           'FontSize',12);
    end
% 
    base1=['perMEAN_' base_name]
    clear base
    base=[base_aux base1 num2str(NFf-NFi+1,5)]
    set(gcf,'PaperPositionMode','auto');
    print('-depsc','-r300',base);
    %print('-djpeg99',base);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if(filme==1)
        j=0;
        for i=1:6
            v(2)=v(2)+10
            base=[base_aux base1 num2str(NFf-NFi+1,5) '_' num2str(j+i-1,5)]
            view(axes1,v);
            set(gcf,'PaperPositionMode','auto');
            print('-depsc','-r300',base);
        end
            for j=1:15
                v(1)=v(1)+10
                base=[base_aux base1 num2str(NFf-NFi+1,5) '_' num2str(j+i-1,5)]
                view(axes1,v);
                set(gcf,'PaperPositionMode','auto');
                print('-depsc','-r300',base);
            end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    clear pcolor
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create figure
figure2 = figure('XVisual',...
    '0xcc (TrueColor, depth 24, RGB mask 0xff0000 0xff00 0x00ff)');

% Create axes
% Create axes
axes2 = axes('Parent',figure2,'XDir','reverse','FontSize',8,...
    'FontName','Times New Roman',...
    'DataAspectRatio',[1 1 3],...%    'CLim',[-2.525 2.562],...
    'CameraViewAngle',12.5);
%view(axes1,[-168.5 -69]);
%view(axes1,[-164.5 -37]);
view(axes2,[150 -45]);
box(axes2,'on');
grid(axes2,'on');
hold(axes2,'all');

slice(X,Y,Z,permmapv,xslice,yslice,zslice), shading flat

% Create xlabel
xlabel('$x$','Interpreter','latex','HorizontalAlignment','right',...
    'FontSize',10,...
    'FontName','Times New Roman');

% Create ylabel
ylabel('$y$','Interpreter','latex','FontSize',10,...
    'FontName','Times New Roman',...
    'HorizontalAlignment','left');

% Create zlabel
zlabel('Time',...%'Interpreter','latex',...
    'FontSize',10,...
    'FontName','Times New Roman',...
    'FontAngle','italic');

% Create colorbar
colorbar('peer',axes2,...
    [0.58 0.35 0.01 0.35],'FontSize',8,...
    'FontName','times');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    posi2=[0.75 0.27 0.045 0.5];
    if(TIPOINPUT==1)
        plot3(vet(:,1),vet(:,2),0.0*vet(:,1),...
            'Parent',axes2,'MarkerFaceColor',[0 0 0],...
            'MarkerEdgeColor',[0 0 0],...
            'MarkerSize',6,...
            'Marker','+','LineStyle','none');
%         %% Create colorbar
%         colorbar1 = colorbar('peer',...
%           axes2,posi2,...
%           'Box','on',...
%           'FontName','times',...
%           'FontSize',12);
    end
    if(TIPOINPUT2==1)
        plot3(vet2(:,1),vet2(:,2),0.0*vet2(:,1),...
            'Parent',axes2,'MarkerFaceColor',[0 0 0],...
            'MarkerEdgeColor',[1 1 1],...
            'MarkerSize',4,...
            'Marker','o','LineStyle','none');
%         %% Create colorbar
%         colorbar1 = colorbar('peer',...
%           axes2,posi2,...
%           'Box','on',...
%           'FontName','times',...
%           'FontSize',12);
    end
% 
    base1=['perVAR_' base_name]
    clear base
    base=[base_aux base1 num2str(NFf-NFi+1,5)]
    set(gcf,'PaperPositionMode','auto');
    print('-depsc','-r300',base);
    %print('-djpeg99',base);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    clear pcolor
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if(stat==1)
        if(grafnorm==1)
            NORMAL(vd,media,std,strtipo)
        else
            LOGNORMAL(vd,media,std,strtipo)
        end
        base1=['histY_' base_name]
        base=[base_aux base1 '_' num2str(II,5)]
        print('-depsc','-r100',base);
    end
    clear vd
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
% close all
