% load_results.m
% charge les fichiers de sortie.
clear
%clf();

escalax=1;
escalay=1;
nb = 3; % engrossar malha do campos de velocidades
AE = 2; % escala dos vetores
disp('----------------')
disp(' LOADING OUTPUT DATA ');

%[FILENAME, PATHNAME] = uigetfile('/home/mrborges/SimuladoresNT/TRACER/sat/*.dat', 'LOAD DATA');
%[FILENAME, PATHNAME] = uigetfile('/home/mrborges/mixed_hybrid/raviart_maicon/ppg-rt-sim-copy/sat/*.dat', 'LOAD DATA');
%[FILENAME, PATHNAME] = uigetfile('/home/mrborges/mixed_hybrid/raviart_maicon/ppg-rt-sim/sat/*.dat', 'LOAD DATA');
[FILENAME, PATHNAME] = uigetfile('../out/ve*.res', 'LOAD DATA');
line_file2=sprintf('%s%s', PATHNAME,FILENAME);
%line_file='../out/sat200x100_hom.res'
%
prompt2={'Diretorio atual figuras/: '};
Answers2=inputdlg(prompt2,'Nome base para os arquivos',1);
base_name = char(Answers2);
base_aux = '../figuras/vel/';
base_name=[base_aux base_name];

time = [];
time1 = [];

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
            xpos = reshape((mattamp(:,1))*escalax,nbx,nby,1);
            ypos = reshape((mattamp(:,2))*escalay,nbx,nby,1);
       end
        vxmat(:,:,nbstep) = reshape(mattamp(:,3),nbx,nby,1);
        vymat(:,:,nbstep) = reshape(mattamp(:,4),nbx,nby,1);
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
MIN=min(min(min((vxmat))));
MAX=max(max(max((vxmat))));
clear vymat vxmat
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('OK, file loaded')
disp(' ')% satmap.m
% cartes de saturation 3D.

disp('----------------')
disp(' 3D saturation movie ');

figure2=figure(...
  'PaperUnits','centimeters',...
  'PaperOrientation','landscape',...
  'PaperSize',[18.92 14.57],...
  'PaperType','<custom>');
%
box on
A = (max(max(xpos)))+dx*0.5;
B = (max(max(ypos)))+dy*0.5;
cc= 0;
dd= 1;
C = cc;
D = 1.1*dd;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nband = 1;
N = length(time1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:N;
%
    quiver(xpos1,ypos1,vx(:,:,i),vy(:,:,i),'LineWidth',1,'AutoScaleFactor',AE,'Color','k')
%
    alpha(0.5);
    shading('interp');
%
    axis([0 A 0 B  C D]);
    daspect([1 1 4*abs(D-C)/A]); 
    set(gca,'ZTick',[cc 0.5*(dd+cc)  dd],'FontName','times','FontSize',12);
%
    zlabel ('S_{w}','FontName','times','FontSize',20);
    set(get(gca,'ZLabel'),'Rotation', 0.0);
    %view(30,15);
    %view(55,23);
    view(0,90);
    num=num2str(i-1);
    base_aux = ['_vel-' num];
    base=[base_name base_aux]
    print('-djpeg90',base)
%    print('-depsc',base)
end

MIN
MAX
clear
disp(' ')
