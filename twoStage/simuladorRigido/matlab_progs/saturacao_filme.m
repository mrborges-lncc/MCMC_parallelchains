% load_results.m
% charge les fichiers de sortie.
clear
%clf();

escalax=1;
escalay=1;

disp('----------------')
disp(' LOADING OUTPUT DATA ');

%[FILENAME, PATHNAME] = uigetfile('/home/mrborges/SimuladoresNT/TRACER/sat/*.dat', 'LOAD DATA');
%[FILENAME, PATHNAME] = uigetfile('/home/mrborges/mixed_hybrid/raviart_maicon/ppg-rt-sim-copy/sat/*.dat', 'LOAD DATA');
%[FILENAME, PATHNAME] = uigetfile('/home/mrborges/mixed_hybrid/raviart_maicon/ppg-rt-sim/sat/*.dat', 'LOAD DATA');
[FILENAME, PATHNAME] = uigetfile('/home/mrborges/SIGER/simulador/out/sa*.res', 'LOAD DATA');
line_file=sprintf('%s%s', PATHNAME,FILENAME);
%line_file='../out/sat200x100_hom.res'
%
prompt2={'Diretorio atual figuras/: '};
Answers2=inputdlg(prompt2,'Nome base para os arquivos',1);
base_name1 = char(Answers2);
base_aux = '../figuras/sat/';
base_name=[base_aux base_name1];

time = [];

fid = fopen(line_file,'r');
test=1;
nbstep = 0;
clear satmat
clear velymat
clear velxmat
MIN=1e+28;
MAX=1e-28;
while test==1
    textline = fgetl(fid);
    if textline == -1
        test=0;
    elseif textline(1:10)=='#TIMESTEP '
        nbstep = nbstep+1;
        time = [time str2num(textline(22:length(textline)))];
        mattamp = fscanf(fid,'%f');
%
        mattamp = reshape(mattamp,3,length(mattamp)/3)';
        if nbstep==1
            dx = 2.0*mattamp(1,1);
            dy = 2.0*mattamp(1,2);
            nbx = int16((max(mattamp(:,1))/dx)+0.5000000001);
            nby = int16((max(mattamp(:,2))/dy)+0.5000000001);
            if(abs(mattamp(1,1)-mattamp(2,1))>1e-7)
                xpos = reshape((mattamp(:,1)),nbx,nby,1);
                ypos = reshape((mattamp(:,2)),nbx,nby,1);
            else
                xpos = reshape((mattamp(:,1)),nby,nbx,1);
                ypos = reshape((mattamp(:,2)),nby,nbx,1);
            end
        end
        if(abs(mattamp(1,1)-mattamp(2,1))>1e-7)
            satmat(:,:,nbstep) = reshape(mattamp(:,3),nbx,nby,1);
        else
            satmat(:,:,nbstep) = reshape(mattamp(:,3),nby,nbx,1);
        end
    end
end
fclose(fid);

disp('OK, file loaded')
disp(' ')% satmap.m
% cartes de saturation 3D.

disp('----------------')
disp(' 3D saturation movie ');
MIN=min(min(min((satmat))));
MAX=max(max(max((satmat))));

figure2=figure(...
  'PaperUnits','centimeters',...
  'PaperOrientation','landscape',...
  'PaperSize',[18.92 14.57],...
  'PaperType','<custom>');
%
M=moviein(length(time));

box on
A = (max(max(xpos)))+dx*0.5;
B = (max(max(ypos)))+dy*0.5;
cc= min(min(min(satmat(:,:,:))));
dd= max(max(max(satmat(:,:,:))));
C = -0.01+cc;
D = 1.1*dd;

for i=1:length(time);
    
    surf(xpos,ypos,satmat(:,:,i));
    hold on
    contour3(xpos,ypos,satmat(:,:,i),20);
    hold off
    shading('interp');
    alpha(0.5);
%
    axis([0 A 0 B  C D]);
    daspect([1 1 4*abs(D-C)/A]); 
    set(gca,'ZTick',[cc 0.5*(dd+cc)  dd],'FontName','times','FontSize',12);
%
    zlabel ('S_{w}','FontName','times','FontSize',20);
    set(get(gca,'ZLabel'),'Rotation', 0.0);
    %view(30,15);
    view(55,23);
    num=num2str(i-1);
    base_aux = ['_sat-' num];
    base=[base_name base_aux]
    print('-djpeg90',base)
%    print('-depsc',base)
    M(i)=getframe;
end

 movie(M,1,1);

test=input(' play again (1/0) ');
while(test==1);
    movie(M,1,3);
    test=input(' play again (1/0) ? ');
end
test = input('store it as AVI (1/0) ? ');
if test==1
%    nomfic = input('file name : ','s');
    nomfic = base_name1;
    nomfic = ['../figuras/movies/' nomfic];
    movie2avi(M,nomfic,'FPS',3);
    movie2avi(M,nomfic,'FPS',4,'quality',100);
end
MIN
MAX
clear
disp(' ')
