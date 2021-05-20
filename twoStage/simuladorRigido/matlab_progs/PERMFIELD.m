%permfield.m
% -----------
close all
escala=1;
clear
prompt2={'Diretorio atual figuras/: '};
Answers2=inputdlg(prompt2,'Nome base para os arquivos',1);
base_name = char(Answers2);
base_aux = '../../figuras/';
%base_aux = '/home/mrborges/Congressos/cilamce2011/figuras/';
%base_aux = '/home/mrborges/THERMOHIDRO/simulador/figuras/';
%base_aux = '/home/mrborges/Dropbox/MCMC_Adaptive/MCMC_DREAM_MOD/figuras/'
Ms = [10000 5000 2500 1000 500 250 200 150 125 100 75 50 25 10 5 1];
Ms = [200 150 210 180];
%
Ms = [0];
N = length(Ms);
imprK=10; % se igual a 1 gera o grafico de K
stat=10;  % se igual a 1 gera os histogramas das vas.
tipo=1 % 1 se perme 2 se phi
permfigs=1; % se igual a 1 gera os graficos 2D
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
M=5.1785e-14;
S=1.0;
M=3.672e+09;
S=0.5;
M=1.0;
S=1.0;
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


disp('----------------')
disp(' LOADING FIELD ');

%[FILENAME, PATHNAME] = uigetfile('../LU/campos/e*.dat', 'LOAD DATA');
%[FILENAME, PATHNAME] =uigetfile('/home/mrborges/RANDOMFIELDS/trunk/FORTRAN/LABTRANGEO_COND/campos/p*.dat', 'LOAD DATA');
%[FILENAME, PATHNAME] =uigetfile('/home/mrborges/RANDOMFIELDS/trunk/OCTAVE/KL/campos/e*.dat', 'LOAD DATA');
%[FILENAME, PATHNAME] =uigetfile('/home/mrborges/RANDOMFIELDS/trunk/OCTAVE/LU/campos/e*.dat', 'LOAD DATA');
%[FILENAME, PATHNAME] =uigetfile('/home/mrborges/paper_KT/simulador_ALE/fieldfrac/kphi*.dat', 'LOAD DATA');
%[FILENAME, PATHNAME] =uigetfile('/home/mrborges/THERMOHIDRO/simulador/field/*.dat', 'LOAD DATA');
%[FILENAME, PATHNAME] =uigetfile('../../gera_KL/MATLAB/campos/*.dat', 'LOAD DATA');
[FILENAME, PATHNAME] =uigetfile('../exp2D/fields/*.dat', 'LOAD DATA');
line_file=sprintf('%s%s', PATHNAME,FILENAME);
%line_file=/home/mrborges/Mcoll/field_klfrac/Yfrac05_s1_M100_120x120_*.dat
cormap = [0 0 0.5625;0 0 0.625;0 0 0.6875;0 0 0.75;0 0 0.8125;0 0 0.875;0 0 0.9375;0 0 1;0 0.0625 1;0 0.125 1;0 0.1875 1;0 0.25 1;0 0.3125 1;0 0.375 1;0 0.4375 1;0 0.5 1;0 0.5625 1;0 0.625 1;0 0.6875 1;0 0.75 1;0 0.8125 1;0 0.875 1;0 0.9375 1;0 1 1;0.0625 1 0.9375;0.125 1 0.875;0.1875 1 0.8125;0.25 1 0.75;0.3125 1 0.6875;0.375 1 0.625;0.4375 1 0.5625;0.5 1 0.5;0.5625 1 0.4375;0.625 1 0.375;0.6875 1 0.3125;0.75 1 0.25;0.8125 1 0.1875;0.875 1 0.125;0.9375 1 0.0625;1 1 0;1 0.9375 0;1 0.875 0;1 0.8125 0;1 0.75 0;1 0.6875 0;1 0.625 0;1 0.5625 0;1 0.5 0;1 0.4375 0;1 0.375 0;1 0.3125 0;1 0.25 0;1 0.1875 0;1 0.125 0;1 0.0625 0;1 0 0;0.9375 0 0;0.875 0 0;0.8125 0 0;0.75 0 0;0.6875 0 0;0.625 0 0;0.5625 0 0;0.5 0 0];
n=0;
pcond = load('../exp2D/input_cond.in');
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

file_base=line_file(1:str_k)
%
ini=str2num(in);
fim=ini+N-1;
for II=ini:1:fim
    line_file = [file_base num2str(Ms(II+1),5) '.dat']
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
            permmap(j,i)=mattamp(k);
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
    ma= 3.0%max(max(permmap))*1.01
    mi=-3.0%min(min(permmap))*1.01
    if(quad==1)
        posi =[0.05 0.15 0.80 0.80];
        posi2=[0.85 0.15 0.05 0.80];
        paperposi=[0.6345 0.6345 11.7 9.7];
    end
    if(quad==0)
        posi=[0.125 0.15 0.76 0.80];
 %       posi=[0.09766 0.3 0.8073 0.8026];
        posi2=[0.925 0.29 0.025 0.52];
        paperposi=[0.6345 0.6345 13 7.5];
    end
    if(quad==2)
        posi=[0.09766 0.3 0.8073 0.8026];
        posi2=[0.925 0.565 0.045 0.27];
    end

    media=mean(mean(permmap))
    vd=reshape(permmap,nx*ny,1);
    variancia=var(vd)
    std=sqrt(variancia)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(permfigs==1)
    figure1 = figure(...
      'PaperUnits','centimeters',...
      'PaperPosition',paperposi,...
      'PaperSize',[12 10],...
      'PaperType','<custom>');

    axes1 = axes('CLim',[mi ma],'Colormap',cormap,...
      'DataAspectRatio',[1 1 2*(ma-mi)/(x*dx)],...%'DataAspectRatio',[x*dx y*dx 3*(ma-mi)/(x*dx)],...
      'FontName','times',...
      'Position',posi,...
      'FontSize',16,...
      'Parent',figure1);
   xlabel(axes1,'$x_1$','FontName','times','FontSize',20,'Interpreter','latex');
   ylabel(axes1,'$x_2$','FontName','times','FontSize',20,'Interpreter','latex');
%     xlabel(axes1,'$$','FontName','times','FontSize',20,'Interpreter','latex');
%     ylabel(axes1,'$$','FontName','times','FontSize',20,'Interpreter','latex');
%    axis(axes1,[0 x*dx 0 y*dx]);
    axis(axes1,[0 x*dx 0 y*dx mi ma]);
    view(axes1,[0 90]);
    %view(axes1,[35 30]);
    box(axes1,'on');
    hold(axes1,'all');
   
    for j=ny:-1:1
        yi = (ny-j)*dy;
        yf = yi+dy;
        yy = [yi yf; yi yf];
        parfor i=1:nx
            xf = i*dx;
            xi = xf - dx;
            xx = [xi xi; xf xf];
            p = permmap(j,i)+[0 0; 0 0];
            surf(xx,yy,p);
            pause(0.005)
        end
    end
    shading('flat')
    % shading('interp')
    % shading('faceted')
    %% Create colorbar
  colorbar(...
      axes1,'Position',posi2,...
      'Box','on',...
      'FontName','times',...
      'FontSize',14);
  hold on
  plot3(pcond(:,1),pcond(:,2),20+pcond(:,3),'LineStyle','none',...
      'LineWidth',1,'Marker','+','Color',[0 0 0],...
      'MarkerSize',12)
% 
    base1=['permg_' base_name]
    clear base
%     base=[base_aux base1 '_' num2str((II+1),5) '_0']
    base=[base_aux base1 '_' num2str(Ms(II+1),5)]
    %set(gcf,'PaperPositionMode','auto');
    print('-depsc','-r300',base);
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
        print('-depsc','-r100',base);
    end
    clear vd
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(imprK==1)
    permmap=M*exp(S*permmap);
    ma=max(max(permmap))
    mi=min(min(permmap))
    media=mean(mean(permmap))
    vd=reshape(permmap,nx*ny,1);
    variancia=var(vd)
    std=sqrt(variancia)
    maa=media+1.96*std
    mii=media-1.96*std
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
        print('-depsc','-r100',base);
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
      'FontSize',16,...
      'Parent',figure2);
%    xlabel(axes2,'$x\ _1$','FontName','times','FontSize',20,'Interpreter','latex');
%    ylabel(axes2,'$x\ _2$','FontName','times','FontSize',20,'Interpreter','latex');
    xlabel(axes2,'$$','FontName','times','FontSize',20,'Interpreter','latex');
    ylabel(axes2,'$$','FontName','times','FontSize',20,'Interpreter','latex');
    axis(axes2,[0 x*dx 0 y*dx mi ma]);
    view(axes2,[0 90]);
    box(axes2,'on');
    hold(axes2,'all');

    for j=ny:-1:1
        yi = (ny-j)*dx;
        yf = yi+dx;
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
    colorbar1 = colorbar('peer',...
      axes2,posi2,...
      'Box','on',...
      'FontName','times',...
      'FontSize',14);
% 
% 
    base1=['permk_' base_name]
    clear base
    base=[base_aux base1 '_' num2str(II,5)]
    set(gcf,'PaperPositionMode','auto');
    print('-depsc','-r100',base);
    %print('-djpeg99',base);
end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
end
clear
