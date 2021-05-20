% permfield.m
% -----------
escala=1;
clear
prompt2={'Diretorio atual figuras/: '};
Answers2=inputdlg(prompt2,'Nome base para os arquivos',1);
base_name = char(Answers2);
base_aux = '../figuras/';
%base_aux = '/home/mrborges/SIGER/simulador/figuras/';

N=1;
imprK=10; % se igual a 1 gera o grafico de K

disp('----------------')
disp(' LOADING FIELD ');

[FILENAME, PATHNAME] = uigetfile('../../campos/po*.dat', 'LOAD DATA');
%[FILENAME, PATHNAME] = uigetfile('/home/mrborges/RANDOMFIELDS/trunk/FORTRAN/LABTRANGEO_COND/campos/p*.dat', 'LOAD DATA');
line_file=sprintf('%s%s', PATHNAME,FILENAME);
%line_file=/home/mrborges/Mcoll/field_klfrac/Yfrac05_s1_M100_120x120_*.dat
n=0;
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


M=1.0;
S=1.0;

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
        quad=0;
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
    figure1 = figure(...
      'PaperUnits','centimeters',...
      'PaperPosition',[0.6345 0.6345 11.7 9.7],...
      'PaperSize',[12 10],...
      'PaperType','<custom>');
  ma=max(max(permmap))
  mi=min(min(permmap))
%      ma=4;
%      mi=-4;
    if(quad==1)
        posi=[0.02 0.1 0.8073 0.8026];
        posi2=[0.75 0.1 0.045 0.8];
    else
        posi=[0.07 0.3 0.8073 0.8026];
        posi2=[0.905 0.49 0.045 0.425];
    end

    axes1 = axes('CLim',[mi ma],...
      'DataAspectRatio',[1 1 2*(ma-mi)/(x*dx)],...%'DataAspectRatio',[x*dx y*dx 3*(ma-mi)/(x*dx)],...
      'FontName','times',...
      'Position',posi,...
      'FontSize',16,...
      'Parent',figure1);
%    xlabel(axes1,'$x\ _1$','FontName','times','FontSize',20,'Interpreter','latex');
%    ylabel(axes1,'$x\ _2$','FontName','times','FontSize',20,'Interpreter','latex');
    xlabel(axes1,'$$','FontName','times','FontSize',20,'Interpreter','latex');
    ylabel(axes1,'$$','FontName','times','FontSize',20,'Interpreter','latex');
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
        for i=1:nx
            xf = i*dx;
            xi = xf - dx;
            xx = [xi xi; xf xf];
            p = permmap(j,i)+[0 0; 0 0];
            surf(xx,yy,p);
        end
    end
    shading('flat')
    % shading('interp')
    % shading('faceted')
    %% Create colorbar
    colorbar1 = colorbar('peer',...
      axes1,posi2,...
      'Box','on',...
      'FontName','times',...
      'FontSize',14);
% 
    base1=['permg_' base_name]
    clear base
    base=[base_aux base1 '_' num2str(II,5)]
    set(gcf,'PaperPositionMode','auto');
    print('-depsc','-r300',base);
    %print('-djpeg99',base);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    clear pcolor
if(imprK==1)
    permmap=M*exp(S*permmap);
    figure2 = figure(...
      'PaperUnits','centimeters',...
      'PaperPosition',[0.6345 0.6345 11.7 9.7],...
      'PaperSize',[12 10],...
      'PaperType','<custom>');
    ma=max(max(permmap))
    mi=min(min(permmap))
    media=mean(mean(permmap))
    std=sqrt(mean(var(permmap)))
    maa=media+std;
    mii=maa-2*std;
    maa=4
    mii=0.4
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
    print('-depsc','-r300',base);
    %print('-djpeg99',base);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
end
clear