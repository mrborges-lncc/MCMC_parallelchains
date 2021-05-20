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
[FILENAME, PATHNAME] = uigetfile('../out/pre*.res', 'LOAD DATA');
line_file=sprintf('%s%s', PATHNAME,FILENAME);
%line_file='../out/pres200x100_hom.res'
prompt2={'Diretorio atual figuras/: '};
Answers2=inputdlg(prompt2,'Nome base para os arquivos',1);
base_name = char(Answers2);
base_aux = '../figuras/pres/';
base_name=[base_aux base_name];


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
%    elseif textline(1:10)=='SATURATION'
        nbstep = nbstep+1;
        time = [time str2num(textline(22:length(textline)))];
        mattamp = fscanf(fid,'%f');
	%nbstep
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
box on
A = (max(max(xpos)))+dx;
B = (max(max(ypos)))+dy;
cc= min(min(min(satmat(:,:,:))));
dd= max(max(max(satmat(:,:,:))));
C = -0.01+cc;
D = 1.1*dd;

for i=1:length(time);
    
    surf(xpos,ypos,satmat(:,:,i));
%    mesh(xpos,ypos,satmat(:,:,i));
    hold on
%    contour3(xpos,ypos,satmat(:,:,i),[0.21 0.268 0.326 0.4 0.5 0.6 0.7 0.8]);
    contour3(xpos,ypos,satmat(:,:,i),20);
%    contour(xpos,ypos,satmat(:,:,i),20);
    hold off
    shading('interp');
    alpha(0.75);
%
% Create axes
    axis([0 A 0 B  C D]);
    daspect([1 1 4*abs(D-C)/A]); 
    set(gca,'ZTick',[cc 0.5*(dd+cc)  dd],'FontName','times','FontSize',12);
    zlabel ('p','FontName','times','FontSize',20);
    set(get(gca,'ZLabel'),'Rotation', 0.0);
    view(55,23);
% %%Create textbox
    num=num2str(i-1);
    base_aux = ['_pre-' num]
    base=[base_name base_aux]
    print('-djpeg90',base)
end

MIN
MAX
clear
disp(' ')
