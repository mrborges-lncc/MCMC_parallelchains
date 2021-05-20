% load_results.m
% charge les fichiers de sortie.
clear
%close all
%clf();
CORTE=1;
escalax=1;
escalay=1;
corte_x=[0 30 70 100]'
% corte_x=[0 20 80 100]'
corte_y=[0 0]'
L=600;
ncov=30;

disp('----------------')
disp(' LOADING OUTPUT DATA ');

%[FILENAME, PATHNAME] = uigetfile('/home/mrborges/SimuladoresNT/TRACER/sat/*.dat', 'LOAD DATA');
%[FILENAME, PATHNAME] = uigetfile('/home/mrborges/mixed_hybrid/raviart_maicon/ppg-rt-sim-copy/sat/*.dat', 'LOAD DATA');
%[FILENAME, PATHNAME] = uigetfile('/home/mrborges/mixed_hybrid/raviart_maicon/ppg-rt-sim/sat/*.dat', 'LOAD DATA');
[FILENAME, PATHNAME] = uigetfile('../simul_comp/exp/out/p*.res', 'LOAD DATA');
%[FILENAME, PATHNAME] = uigetfile('/home/mrborges/SIGER/simuladorCO2/out/sa*.res', 'LOAD DATA');
line_file=sprintf('%s%s', PATHNAME,FILENAME);
%line_file='../out/sat200x100_hom.res'
%
% prompt2={'Diretorio atual figuras/: '};
% Answers2=inputdlg(prompt2,'Nome base para os arquivos',1);
% base_name = char(Answers2);
base_aux = '../figuras/';
% base_name=[base_aux base_name];

cl1=5e8;      % minimo valor no colormap
cl2=1.5e7;      % maximo valor no colormap
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
nt=size(time,2)-1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% CORTE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(CORTE==1)
    cortex=int16(corte_x/dx);
    cortey=int16(corte_y/dy);
    numx=cortex(2)-cortex(1)+(cortex(4)-cortex(3));
    newsat=zeros(numx,cortey(1)+(nby-cortey(2)),nt);
    x=zeros(numx,1);
    for j=1:cortex(2)-cortex(1)
        x(j)=cortex(1)*dx+dx/2+(j-1)*dx;
    end
    for i=1:cortex(4)-cortex(3)
        x(j+i)=(cortex(3))*dx+dx/2+(i-1)*dx;
    end
    y=zeros(cortey(1)+(nby-cortey(2)),1);
    for i=1:cortey(1)
        y(i)=(i-1)*dy+dy/2;
    end
    for i=1:nby-cortey(2)
        y(i+cortey(1))=(cortey(2)+(i-1))*dy+dy/2;
    end
    [Y,X]=meshgrid(y,x);
    for k=1:nt
        for i=1:numx
            for j=1:cortey(1)+(nby-cortey(2))
                xp=X(i,j);
                yp=Y(i,j);
                for jj=1:nby
                    for ii=1:nbx
                        if((abs(xp-xpos(ii,jj))<1e-5)&&...
                                (abs(yp-ypos(ii,jj)))<1e-5)
                            newsat(i,j,k)=satmat(ii,jj,k);
                        end
                    end
                end
            end
        end
    end
    nbx=size(newsat,1);
    nby=size(newsat,2);
    satmat=newsat;
    clear newsat
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% COVARIANCIA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
NN=4;
ntN=int32(nt/NN)
cov =zeros(ntN,1);
t   =zeros(ntN,1);
cont=zeros(ntN,1);
med =mean(mean(mean(satmat)));
for k=0:ntN-1
    k
    for i=1:nbx
        for j=1:nby
            for kk=1:nt-k
                medkk =mean(mean(satmat(:,:,kk)));
                medkkk=mean(mean(satmat(:,:,kk+k)));
                cov(k+1)=cov(k+1)+(satmat(i,j,kk)-medkk)*...
                    (satmat(i,j,kk+k)-medkkk);
                cont(k+1)=cont(k+1)+1;
                t(k+1) = time(kk+k)-time(kk);
            end
        end
    end
end
cov=cov./cont;
cov=cov/cov(1,1);
%
%if(ncov==3)
    L=figregress(t,log(cov))
    co=exp(-0.5*(t.^2)/L^2);
    figCVP(t,[cov co]);
%else
    L1=figregress1(t,log(cov))
    co=exp(-(t.^1)/L1^1);
    figCVP(t,[cov co]);
%end
disp('----------------')
disp(' 3D saturation movie ');
MIN=min(min(min((satmat))));
MAX=max(max(max((satmat))));
%
box on
A = 100;%(max(max(X)))+dx*0.5;
B = (max(max(Y)))+dy*0.5;
aa = 0.0;
bb = 0.0;
cc= min(min(min(satmat(:,:,:))));
dd= max(max(max(satmat(:,:,:))));
cl1=cc;
cl2=dd;
% cc= cl1;
% dd= cl2;
C = 0.9*cc;
D = 1.2*dd;
 
for i=1:1%length(time);
    t=time(i)
    num=num2str(i-1);
    if(contorno==0)
        sat_figure(X,Y,satmat(:,:,i),A,B,aa,bb,cc,dd,C,D,t,cl1,cl2);
        base_aux = ['_pre-' num];
        base=[base_name base_aux];
%         print('-djpeg90',base)
       print('-depsc','-r100',base)
    else
        sat_figure_contour(X,Y,satmat(:,:,i),A,B,aa,bb,cc,dd,C,D,t,cl1,cl2);
        base_aux = ['_prec-' num];
        base=[base_name base_aux];
        print('-djpeg90',base)
%        print('-depsc','-r300',base)
    end
end
%close all
MIN
MAX
%clear
disp(' ')
% close all