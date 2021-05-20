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
[FILENAME, PATHNAME] = uigetfile('../out/b*.res', 'LOAD DATA ONE');
%[FILENAME, PATHNAME] = uigetfile('/home/mrborges/SIGER/simuladorCO2/out/sa*.res', 'LOAD DATA');
line_file=sprintf('%s%s', PATHNAME,FILENAME);
[FILENAME, PATHNAME] = uigetfile('../out/b*.res', 'LOAD DATA TWO');
%[FILENAME, PATHNAME] = uigetfile('/home/mrborges/SIGER/simuladorCO2/out/sa*.res', 'LOAD DATA');
line_file1=sprintf('%s%s', PATHNAME,FILENAME);
%line_file='../out/sat200x100_hom.res'
%
% prompt2={'Diretorio atual figuras/: '};
% Answers2=inputdlg(prompt2,'Nome base para os arquivos',1);
% base_name = char(Answers2);
base_aux = '../figuras/mass/';
% base_name=[base_aux base_name];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=length(FILENAME):-1:1
    a=FILENAME(i);
    if(a=='.')
        LF=i-1;
    end
end
base_name = [base_aux FILENAME(5:LF)];

time1 = [];

fid = fopen(line_file1,'r');
test=1;
nbstep = 0;
clear satmat1
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
        time1 = [time1 str2num(textline(22:length(textline)))];
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
            satmat1(:,:,nbstep) = reshape(mattamp(:,3),nbx,nby,1);
        else
            satmat1(:,:,nbstep) = reshape(mattamp(:,3),nby,nbx,1);
        end
    end
end
fclose(fid);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('OK, file loaded')
disp(' ')% satmap.m
% cartes de saturation 3D.

disp('----------------')
MIN=min(min(min((satmat))))
MAX=max(max(max((satmat))))
MIN1=min(min(min((satmat1))))
MAX1=max(max(max((satmat1))))

%
MB=time*0.0;
MB1=time1*0.0;
for i=2:length(time);
    MB(i)=max(max(abs((satmat(:,:,i)))));
end
for i=2:length(time1);%     else
%         mass_figure_contour(xpos,ypos,satmat(:,:,i),A,B,aa,bb,cc,dd,C,D,t,cl1,cl2);
%         base_aux = ['_massMAXc-' num];
%         base=[base_name base_aux];
% %        print('-djpeg90',base)
%         print('-depsc','-r300',base)
%     end
% end
% close all

    MB1(i)=max(max(abs((satmat1(:,:,i)))));
end
% pl1=line(time(2:end),MB(2:end),'Color','r');
% ax1=gca;
% set(ax1,'XColor','r','YColor','r');
% ax2 = axes('Position',get(ax1,'Position'),...
%            'XAxisLocation','top',...
%            'YAxisLocation','right',...
%            'Color','none',...
%            'XColor','k','YColor','k');
% pl2=line(time1(2:end),MB1(2:end),'Color','b','Parent',ax2);
 maxerrorfig(time(2:end),MB(2:end),time1(2:end),MB1(2:end))     

base=[base_aux 'MAX' FILENAME(5:LF)];
%        print('-djpeg90',base)
print('-depsc','-r300',base)
MIN
MAX
clear
disp(' ')
