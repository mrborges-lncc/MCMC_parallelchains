%permfield.m
% -----------
escala=1;
%close all;
clear;
%prompt2={'Diretorio atual figuras/: '};
%Answers2=inputdlg(prompt2,'Nome base para os arquivos',1);
%base_name = char(Answers2);
%base_aux = '../figuras/';
%base_aux = '/home/mrborges/Congressos/cilamce2011/figuras/';
%base_aux = '/home/mrborges/THERMOHIDRO/simulador/figuras/';
base_aux = '../figuras/'
N=1;
imprK=10; % se igual a 1 gera o grafico de K
stat=10;  % se igual a 1 gera os histogramas das vas.
tipo=1 % 1 se perm e 2 se phi
permfigs=1; % se igual a 1 gera os graficos 2D
TIPOINPUT = 1; % ==1 entrada por arquivo
TIPOINPUT2 = 1; % ==1 entrada por arquivo
mi=-3.5;
ma=3.5;
controle=0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
M=7.5e-14;
S=0.73;
%S=0.225;
%M=0.215;
M=3.475e09
S=0.25
M=2.87e08
S=0.2
M= 0.197546
S= 0.157139
%
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
inp = load('../gera_KL/FORTRAN_RW/in/input_cond.dat');
%
if(TIPOINPUT2 == 1)
    name = '../simuladorRigido/exp/input_cond.in'; 
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

%[FILENAME, PATHNAME] = uigetfile('../LU/campos/e*.dat', 'LOAD DATA');
%[FILENAME, PATHNAME] =uigetfile('/home/mrborges/RANDOMFIELDS/trunk/FORTRAN/LABTRANGEO_COND/campos/p*.dat', 'LOAD DATA');
%[FILENAME, PATHNAME] =uigetfile('/media/Dados/RANDOMFIELDS/trunk/OCTAVE/KL/campos/*.dat', 'LOAD DATA');
%[FILENAME, PATHNAME] =uigetfile('/home/mrborges/RANDOMFIELDS/trunk/OCTAVE/LU/campos/e*.dat', 'LOAD DATA');
%[FILENAME, PATHNAME] =uigetfile('/home/mrborges/paper_KT/simulador_ALE/fieldfrac/kphi*.dat', 'LOAD DATA');
%[FILENAME, PATHNAME] =uigetfile('/home/mrborges/THERMOHIDRO/simulador/field/*.dat', 'LOAD DATA');
%[FILENAME, PATHNAME] =uigetfile('../SIMULADOR_VISCOELASTICO/exp/fields/*.dat', 'LOAD DATA');
%[FILENAME, PATHNAME] =uigetfile('../twoStage/SIMULADOR_VISCOELASTICO/exp/fields/p*.dat', 'LOAD DATA');
[FILENAME, PATHNAME] =uigetfile('../simuladorRigido/exp000/fields/*.dat', 'LOAD DATA');
% [FILENAME, PATHNAME] =uigetfile('../twoStage/select_fields/*.dat', 'LOAD DATA');
%[FILENAME, PATHNAME] =uigetfile('~/TwoWay_lastversionPAPER/trunk/SIMULADOR_VISCOELASTICO/expSISMELAST/fields/*.dat', 'LOAD DATA');
line_file=sprintf('%s%s', PATHNAME,FILENAME);
%line_file='../simuladorBTMM/exp01/fields/ref_0.dat'
n=0;
%
for i=size(line_file,2):-1:1
   a=line_file(i);
   if(a=='.')
       in=line_file(i-1);
       for j=i-2:-1:1
           a=line_file(j);
           if(a=='_')
               str_k=j;
               i=0;
               break
           else
               in=[a in]
           end
       end
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
ini=str2num(in);
fim=ini+N-1;
for II=ini:1:fim
    line_file = [file_base num2str(II,5) '.dat']
    fid = fopen(line_file,'r');
    mattamp = fscanf(fid,'%f');

    disp('file loaded.')
    fclose(fid);
    inf = mattamp(1:4);
    Lx = inf(1);
    Ly = inf(2);
    nx = inf(3);
    ny = inf(4);
    dx = Lx/nx;
    dy = Ly/ny;
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
    permmap=[];
    k=0;
    for j=ny:-1:1
        k=k+1;
        if(mattamp(k)~=ny-j)
            disp('erro1')
            break
        end
        for i=1:nx
            k=k+1;
            permmap(j,i)=S*mattamp(k);
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
        ma=max(max(permmap))*1.01
        mi=min(min(permmap))*1.01
    end
    if(quad==1)
        posi=[0.02 0.1 0.8073 0.8026];
        posi2=[0.75 0.1 0.045 0.8];
    end
    if(quad==0)
        posi=[0.0976 0.3 0.76 0.8026];
 %       posi=[0.09766 0.3 0.8073 0.8026];
        posi2=[0.87 0.49 0.045 0.425];
    end
    if(quad==2)
        posi=[0.09766 0.3 0.8073 0.8026];
        posi2=[0.925 0.565 0.045 0.27];
        posi2=[0.93 0.64 0.01 0.115];
    end

    media=mean(mean(permmap))
    vd=reshape(permmap,nx*ny,1);
    variancia=var(vd)
    std=sqrt(variancia)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(permfigs==1)
    figure1 = figure(...
      'PaperUnits','centimeters',...
      'PaperPosition',[0.6345 0.6345 11.7 5.7],...
      'PaperSize',[12 6],...
      'PaperType','<custom>');

    axes1 = axes('CLim',[mi ma],...
      'DataAspectRatio',[1 1 2*(ma-mi)/(x*dx)],...%'DataAspectRatio',[x*dx y*dx 3*(ma-mi)/(x*dx)],...
      'FontName','times',...
      'Position',posi,...
      'FontSize',12,...
      'Parent',figure1);
%    xlabel(axes1,'$x\ _1$','FontName','times','FontSize',20,'Interpreter','latex');
%    ylabel(axes1,'$x\ _2$','FontName','times','FontSize',20,'Interpreter','latex');
    xlabel(axes1,'$x (m)$','FontName','times','FontSize',16,'Interpreter','latex');
    ylabel(axes1,'$y (m)$','FontName','times','FontSize',16,'Interpreter','latex');
%    axis(axes1,[0 x*dx 0 y*dx]);
    axis(axes1,[0 x*dx 0 y*dy mi ma]);
    view(axes1,[0 90]);
    %view(axes1,[35 30]);
    box(axes1,'on');
    hold(axes1,'all');
   
    for j=ny:-1:1
        yi = (ny-j)*dy;
        yf = yi+dy;
        yy = [yi yf; yi yf];
        for i=1:nx
            xf = i*dx;
            xi = xf - dx;
            xx = [xi xi; xf xf];
            p = permmap(j,i)+[0 0; 0 0];
            surf(xx,yy,p);
        end
    end
    shading('flat')
    %shading('interp')
    % shading('faceted')
    if(TIPOINPUT==1)
        plot3(vet(:,1),vet(:,2),ma+0.0*vet(:,1),...
            'Parent',axes1,'MarkerFaceColor',[0 0 0],...
            'MarkerEdgeColor',[0 0 0],...
            'MarkerSize',14,...
            'Marker','+','LineStyle','none');
    end
    if(TIPOINPUT2==1)
        plot3(vet2(:,1),vet2(:,2),ma+0.0*vet2(:,1),...
            'Parent',axes1,'MarkerFaceColor',[0 0 0],...
            'MarkerEdgeColor',[1 1 1],...
            'MarkerSize',6,...
            'Marker','o','LineStyle','none');
    end
        %% Create colorbar
        colorbar1 = colorbar(...
          'Box','on',...
          'FontName','times',...
          'FontSize',10);
% 
    base1=['permg_' base_name]
    clear base
    base=[base_aux base1 num2str(II,5)]
    set(gcf,'PaperPositionMode','auto');
    print('-depsc','-r600',base);
    %print('-djpeg99',base);
end
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
        print('-depsc','-r600',base);
    end
    clear vd
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(imprK==1)
    disp('KKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKK')
    permmap=M*exp(permmap);
    ma=max(max(permmap))
    mi=min(min(permmap))
    media=mean(mean(permmap))
    vd=reshape(permmap,nx*ny,1);
    variancia=var(vd)
    std=sqrt(variancia)
    maa=media+1.96*std;
    mii=media-1.96*std;
    if(mii<0)
        mii=mi
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if(stat==1)
        if(grafnorm==0)
            NORMAL(vd,media,std,strtipo)
        else
            LOGNORMAL(vd,media,std,strtipo)
        end
        base1=['histK_' base_name]
        base=[base_aux base1 '_' num2str(II,5)]
        print('-depsc','-r600',base);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(permfigs==1)
    figure2 = figure(...
      'PaperUnits','centimeters',...
      'PaperPosition',[0.6345 0.6345 11.7 9.7],...
      'PaperSize',[12 10],...
      'PaperType','<custom>');
    axes2 = axes('CLim',[mii maa],...
      'DataAspectRatio',[1 1 2*(ma-mi)/(x*dx)],...%'DataAspectRatio',[x*dx y*dx 3*(ma-mi)/(x*dx)],...
      'FontName','times',...
      'Position',posi,...
      'FontSize',12,...
      'Parent',figure2);
%    xlabel(axes2,'$x\ _1$','FontName','times','FontSize',20,'Interpreter','latex');
%    ylabel(axes2,'$x\ _2$','FontName','times','FontSize',20,'Interpreter','latex');
    xlabel(axes2,'$x$','FontName','times','FontSize',16,'Interpreter','latex');
    ylabel(axes2,'$y$','FontName','times','FontSize',16,'Interpreter','latex');
    axis(axes2,[0 x*dx 0 y*dy mi ma]);
    view(axes2,[0 90]);
    box(axes2,'on');
    hold(axes2,'all');

    for j=ny:-1:1
        yi = (ny-j)*dy;
        yf = yi+dy;
        yy = [yi yf; yi yf];
        for i=1:nx
            xf = i*dx;
            xi = xf - dx;
            xx = [xi xi; xf xf];
            p = permmap(j,i)+[0 0; 0 0];
            surf(xx,yy,p);
        end
    end
        ma=max(max(permmap))*1.01

    shading('flat')
%     colorbar1 = colorbar('peer',...
%       axes2,posi2,...
%       'Box','on',...
%       'FontName','times',...
%       'FontSize',14);
% % 
%     plot3(vet(:,1),vet(:,2),ma+0.0*vet(:,1),...
%         'Parent',axes2,'MarkerFaceColor','none',...
%         'MarkerEdgeColor',[0 0 0],...
%         'MarkerSize',12,...
%         'Marker','+','LineStyle','none');
     if(TIPOINPUT==1)
        plot3(vet(:,1),vet(:,2),ma+0.0*vet(:,1),...
            'Parent',axes2,'MarkerFaceColor',[0 0 0],...
            'MarkerEdgeColor',[0 0 0],...
            'MarkerSize',14,...
            'Marker','+','LineStyle','none');
     end
    if(TIPOINPUT2==1)
        plot3(vet2(:,1),vet2(:,2),ma+0.0*vet2(:,1),...
            'Parent',axes2,'MarkerFaceColor',[0 0 0],...
            'MarkerEdgeColor',[1 1 1],...
            'MarkerSize',6,...
            'Marker','o','LineStyle','none');
    end
        %% Create colorbar
        colorbar2 = colorbar('peer',...
          axes2,posi2,...
          'Box','on',...
          'FontName','times',...
          'FontSize',10);

    base1=['permk_' base_name]
    clear base
    base=[base_aux base1 num2str(II,5)]
    set(gcf,'PaperPositionMode','auto');
    print('-depsc','-r600',base);
    %print('-djpeg99',base);
end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
end
clear
%close all
