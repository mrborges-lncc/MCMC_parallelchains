% load_results.m
% charge les fichiers de sortie.
clear
%clf();

escalax=1;
escalay=1;
nb = 5; % engrossar malha do campos de velocidades
AE = 16/nb; % escala dos vetores
disp('----------------')
disp(' LOADING OUTPUT DATA ');

%[FILENAME, PATHNAME] = uigetfile('/home/mrborges/SimuladoresNT/TRACER/sat/*.dat', 'LOAD DATA');
%[FILENAME, PATHNAME] = uigetfile('/home/mrborges/mixed_hybrid/raviart_maicon/ppg-rt-sim-copy/sat/*.dat', 'LOAD DATA');
%[FILENAME, PATHNAME] = uigetfile('/home/mrborges/mixed_hybrid/raviart_maicon/ppg-rt-sim/sat/*.dat', 'LOAD DATA');
[FILENAME, PATHNAME] = uigetfile('../out/sa*.res', 'LOAD DATA');
line_file=sprintf('%s%s', PATHNAME,FILENAME);
[FILENAME, PATHNAME] = uigetfile('../out/ve*.res', 'LOAD DATA');
line_file2=sprintf('%s%s', PATHNAME,FILENAME);
%line_file='../out/sat200x100_hom.res'
%
prompt2={'Diretorio atual figuras/: '};
Answers2=inputdlg(prompt2,'Nome base para os arquivos',1);
base_name = char(Answers2);
base_aux = '../figuras/sat/';
base_name=[base_aux base_name];

time = [];
time1 = [];

fid = fopen(line_file,'r');
test=1;
nbstep = 0;
clear satmat
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear mattamp
fid = fopen(line_file2,'r');
test=1;
nbstep = 0;
%
while test==1
    textline = fgetl(fid);
    if textline == -1
        test=0;
    elseif textline(1:10)=='#TIMESTEP '
        nbstep = nbstep+1;
        time1 = [time1 str2num(textline(22:length(textline)))];
        mattamp = fscanf(fid,'%f');
%
        mattamp = reshape(mattamp,4,length(mattamp)/4)';
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
            vxmat(:,:,nbstep) = reshape(mattamp(:,3),nbx,nby,1);
            vymat(:,:,nbstep) = reshape(mattamp(:,4),nbx,nby,1);
        else
            vxmat(:,:,nbstep) = reshape(mattamp(:,3),nby,nbx,1);
            vymat(:,:,nbstep) = reshape(mattamp(:,4),nby,nbx,1);
        end
    elseif textline(1:10)=='VELOCITY -'
        mattamp = fscanf(fid,'%f');
        mattamp = reshape(mattamp,4,length(mattamp)/4)';
    end
end
fclose(fid);
clear mattamp
vx = vxmat(1:nb:end,1:nb:end,:);
vy = vymat(1:nb:end,1:nb:end,:);
xpos1 = xpos(1:nb:end,1:nb:end,:);
ypos1 = ypos(1:nb:end,1:nb:end,:);
clear vymat vxmat
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
%M=moviein(length(time));

box on
A = (max(max(xpos)))+dx*0.5;
B = (max(max(ypos)))+dy*0.5;
cc= min(min(min(satmat(:,:,:))));
dd= max(max(max(satmat(:,:,:))));
C = 0;-0.1+cc;
D = 1.1*dd;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nband = 1;
N = min(length(time),length(time1));
if(length(time1)==2);
    nband = 0;
    N = length(time);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:N;
%    
    if(nband==1)
        surf(xpos,ypos,satmat(:,:,i));
        hold on
        quiver(xpos1,ypos1,vx(:,:,i),vy(:,:,i),'LineWidth',1,...
            'AutoScaleFactor',AE,'Color','k');
        hold off
    else
        surf(xpos,ypos,satmat(:,:,i));
        hold on
        quiver(xpos1,ypos1,vx(:,:,2),vy(:,:,2),'LineWidth',1,...
            'AutoScaleFactor',AE,'Color','k');
        hold off
    end
%
    alpha(0.5);
    shading('interp');
    axis([0 A 0 B  C D]);
    daspect([1 1 4*abs(D-C)/A]); 
    set(gca,'ZTick',[cc 0.5*(dd+cc)  dd],'FontName','times','FontSize',12);
%
    zlabel ('S_{w}','FontName','times','FontSize',20);
    set(get(gca,'ZLabel'),'Rotation', 0.0);
    view(0,90);
%
    num=num2str(i-1);
    base_aux = ['_satvel-' num];
    base=[base_name base_aux]
    print('-djpeg90',base)
%    print('-depsc',base)
end

MIN
MAX
clear
disp(' ')
