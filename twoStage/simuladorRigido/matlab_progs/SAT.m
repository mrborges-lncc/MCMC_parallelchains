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
[FILENAME, PATHNAME] = uigetfile('../experimentos/exp2D/out/s*.res', 'LOAD DATA');
%[FILENAME, PATHNAME] = uigetfile('/home/mrborges/SIGER/simuladorCO2/out/sa*.res', 'LOAD DATA');
line_file=sprintf('%s%s', PATHNAME,FILENAME);
%line_file='../out/sat200x100_hom.res'
%
% prompt2={'Diretorio atual figuras/: '};
% Answers2=inputdlg(prompt2,'Nome base para os arquivos',1);
% base_name = char(Answers2);
base_aux = '../figuras/';
% base_name=[base_aux base_name];

cl1=0.0;      % minimo valor no colormap
cl2=1.0;      % maximo valor no colormap
contorno = 0; % se diferente de zero graf. contorno 
for i=length(FILENAME):-1:1
    a=FILENAME(i);
    if(a=='.')
        LF=i-1;
    end
end
base_name = [base_aux FILENAME(5:LF)];

time = [];

fid = fopen(line_file,'r');
test=1;
nbstep = 0;
clear satmat
clear velymat
clear velxmat
MIN=-1e+28;
MAX=1e+28;
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

cl11=0.01;
TOL=1e-3;
%%AJUSTE PARA A FIGURA
% for k=1:length(time)
%     for i=1:nbx
%         for j=1:nby
%             if(abs(satmat(i,j,k)-cl11)<TOL)
%                 satmat(i,j,k)=cl11;
%             end
%         end
%     end
% end

% cartes de saturation 3D.

disp('----------------')
disp(' 3D saturation movie ');
MIN=min(min(min((satmat))))
MAX=max(max(max((satmat))))

%
box on
A = (max(max(xpos)))+dx*0.5;
B = (max(max(ypos)))+dy*0.5;
aa = 0.0;
bb = 0.0;
% cc= min(min(min(satmat(:,:,:))));
% dd= max(max(max(satmat(:,:,:))));
cc= cl1;
dd= cl2;
C = 0.9*cc;
D = 1.2*dd;
C = cc;
D = dd;

for i=1:length(time);
    t=time(i)
    num=num2str(i-1);
    if(contorno==0)
        sat_figure(xpos,ypos,satmat(:,:,i),A,B,aa,bb,cc,dd,C,D,t,cl1,cl2);
        base_aux = ['_sat-' num];
        base=[base_name base_aux];
        print('-zbuffer','-djpeg90',base)
%        print('-zbuffer','-depsc','-r100',base)
%        print('-depsc','-r100',base)
    else
        sat_figure_contour(xpos,ypos,satmat(:,:,i),A,B,aa,bb,cc,dd,C,D,t,cl1,cl2);
        base_aux = ['_satc-' num];
        base=[base_name base_aux];
%        print('-djpeg90',base)
        print('-depsc','-r100',base)
    end
end
poro=0.2;
massaTotal = sum(sum(satmat(:,:,end)))*poro*dx*dy
close all
MIN
MAX
clear
disp(' ')
