%permfield.m
% -----------
escala=1;
clear;
%prompt2={'Diretorio atual figuras/: '};
%Answers2=inputdlg(prompt2,'Nome base para os arquivos',1);
%base_name = char(Answers2);
%base_aux = '../figuras/';
%base_aux = '/home/mrborges/Congressos/cilamce2011/figuras/';
%base_aux = '/home/mrborges/THERMOHIDRO/simulador/figuras/';
base_aux = '../figuras/'
nz=100;
NFi=0;
NFf=99;
Lz=1000;
Lx=100;
Ly=100;
nx=50;
ny=50;
nz=100;
dx=Lx/nx;
dy=Ly/ny;
dz=Lz/nz;
imprK=10; % se igual a 1 gera o grafico de K
stat=10;  % se igual a 1 gera os histogramas das vas.
tipo=1 % 1 se perm e 2 se phi
permfigs=1; % se igual a 1 gera os graficos 2D
TIPOINPUT = 10; % ==1 entrada por arquivo
TIPOINPUT2 = 10; % ==1 entrada por arquivo
filme=10;
mi=-3;
ma=3;
controle=0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
M=4.7887e-14;
S=1.0;
M=1.0;
% S=0.2918;
% M=0.2396;
%S=1.0
%M=1.0
% M=-1.3671
% S=0.0791
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
inp = load('../gera_KL/MATLAB/input_cond.dat');
%
if(TIPOINPUT2 == 1)
    name = '../simul_comp/exp/input_cond.in'
%     name = '../simuladorBTMM/exp/input_cond.in'; 
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

% [FILENAME, PATHNAME] =uigetfile('../forecast/3Dfields/fields/*.dat', 'LOAD DATA');
% [FILENAME, PATHNAME] =uigetfile('../simul_comp/exp/fields/*.dat', 'LOAD DATA');
% [FILENAME, PATHNAME] =uigetfile('../simuladorBTMM/exp/fields/*.dat', 'LOAD DATA');
% [FILENAME, PATHNAME] =uigetfile('../twoStage/select_fields/*.dat', 'LOAD DATA');
% [FILENAME, PATHNAME] =uigetfile('../gera_KL/MATLAB/campos/*.dat', 'LOAD DATA');
% line_file=sprintf('%s%s', PATHNAME,FILENAME);
line_file='../forecast/3Dfields/fields/MC_0_0.dat'
line_file='../forecast/3Dfields/fields/MC_0_0.dat'
line_file='../forecast/FORTRAN_LBD3D/campos/MCMC_0_0.dat'
n=0;
%
% for i=size(line_file,2):-1:1
%    a=line_file(i);
%    if(a=='.')
%        in=line_file(i-1);
%        for j=i-2:-1:1
%            a=line_file(j);
%            if(a=='_')
%                str_k=j;
%                i=0;
%                break
%            else
%                in=[a in]
%            end
%        end
%    end
% end
sgn=0;
for i=size(line_file,2):-1:1
    a=line_file(i);
    if(a=='_')
        sgn=sgn+1;
        line_file(1:i)
    end
    if(sgn==2)
        str_k=i
        break
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
ini=0;
fim=ini+nz-1;
for nf=NFi:NFf
    cont=0;
    for II=ini:1:fim
        line_file = [file_base num2str(nf,5) '_' num2str(II,5) '.dat']
        fid = fopen(line_file,'r');
        mattamp = fscanf(fid,'%f');
        cont=cont+1;
        disp('file loaded.')
        fclose(fid);
        inf = mattamp(1:4);
        Lx = inf(1);
        Ly = inf(2);
        nx = inf(3);
        ny = inf(4);
        dx = Lx/nx;
        dy = Ly/ny;
        if(nf==NFi&&II==ini)
            permmap =zeros(ny,nx,nz);
            permmapv=zeros(ny,nx,nz);
        end
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
        k=0;
        for j=ny:-1:1
            k=k+1;
            if(mattamp(k)~=ny-j)
                disp('erro1')
                break
            end
            for i=1:nx
                k=k+1;
                permmap(j,i,cont) =permmap(j,i,cont)+(S*mattamp(k));
    %            permmapv(j,i,cont)=permmapv(j,i,cont)+(S*mattamp(k))^2;
    %             permmap(j,i,cont)=M*exp(S*mattamp(k));
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
            ma=max(max(permmap))*1.01;
            mi=min(min(permmap))*1.01;
        end
        clear vd mattamp
    end
end
media    =permmap/(NFf-NFi+1);
media = 0.0*media;
%variancia=permmapv/(NFf-NFi+1)-permmap.^2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
covx=zeros(nx/2,1);
covy=zeros(ny/2,1);
covz=zeros(nz/2,1);
contx=zeros(nx/2,1);
conty=zeros(ny/2,1);
contz=zeros(nz/2,1);
for nf=NFi:NFf
    cont=0;
    for II=ini:1:fim
        line_file = [file_base num2str(nf,5) '_' num2str(II,5) '.dat']
        fid = fopen(line_file,'r');
        mattamp = fscanf(fid,'%f');
        cont=cont+1;
        disp('file loaded.')
        fclose(fid);
        inf = mattamp(1:4);
        Lx = inf(1);
        Ly = inf(2);
        nx = inf(3);
        ny = inf(4);
        dx = Lx/nx;
        dy = Ly/ny;
        if(nf==NFi&&II==ini)
            permmap =zeros(ny,nx,nz);
            permmapv=zeros(ny,nx,nz);
        end
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
        k=0;
        for j=ny:-1:1
            k=k+1;
            if(mattamp(k)~=ny-j)
                disp('erro1')
                break
            end
            for i=1:nx
                k=k+1;
                permmap(j,i,cont) =(S*mattamp(k));
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
            ma=max(max(permmap))*1.01;
            mi=min(min(permmap))*1.01;
        end
        clear vd mattamp
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% covariancia na direcao x %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for i=1:nx/2
        d=i-1;
        for z=1:nz
            for k=1:ny
                for j=1:nx-d
                    covx(i)=covx(i)+(permmap(k,j,z)-media(k,j,z))*...
                        (permmap(k,j+d,z)-media(k,j+d,z));
                    contx(i)=contx(i)+1;
                end
            end
        end
    end            
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% covariancia na direcao y %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for i=1:ny/2
        d=i-1;
        for z=1:nz
            for k=1:nx
                for j=1:ny-d
                    covy(i)=covy(i)+(permmap(j,k,z)-media(j,k,z))*...
                        (permmap(j+d,k,z)-media(j+d,k,z));
                    conty(i)=conty(i)+1;
                end
            end
        end
    end            
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% covariancia na direcao z %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for i=1:nz/2
        d=i-1;
        for z=1:ny
            for k=1:nx
                for j=1:nz-d
                    covz(i)=covz(i)+(permmap(z,k,j)-media(z,k,j))*...
                        (permmap(z,k,j+d)-media(z,k,j+d));
                    contz(i)=contz(i)+1;
                end
            end
        end
    end            
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
clear media permmap
covx=covx./contx;
covy=covy./conty;
covz=covz./contz;
x=[0:dx:Lx/2-dx]';
y=[0:dy:Ly/2-dy]';
z=[0:dz:Lz/2-dz]';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
X=[x covx];
Y=[y covy];
Z=[z covz];
home='/Users/mrborges/projMCMC/trunk/forecast/FORTRAN_LBD3D/out/';
nomex= [home 'gx.dat']
nomey= [home 'gy.dat']
nomez= [home 'gz.dat']

save(nomex,'X','-ascii');
save(nomey,'Y','-ascii');
save(nomez,'Z','-ascii');
clear
