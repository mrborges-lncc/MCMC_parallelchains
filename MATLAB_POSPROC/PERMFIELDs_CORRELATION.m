%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FIELD CORRELATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
file_base = '../twoStage/select_fields/'
file_base = '~/fields/campos/';
[FILENAME, PATHNAME] =uigetfile([file_base '*.dat'], 'LOAD DATA');
nc=getNAME(FILENAME);
namet1 = FILENAME(1:nc);
[FILENAME, PATHNAME] =uigetfile([file_base '*.dat'], 'LOAD DATA');
nc=getNAME(FILENAME);
namet2 = FILENAME(1:nc);
%namet1 = 'field_v1_KL20_'
%namet2 = 'field_v2_KL20_'
ini=60;
fim=60;
N =fim-ini+1;
Sk=1.0;
Mk=4.8e-14; 
Sp=0.2864;
Mp=0.2419;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('----------------')
disp(' LOADING FIELD ');
%
base_name = [file_base namet1];
%
num=0;
fim=ini+N-1;
for II=ini:1:fim
    line_file = [base_name num2str(II,5) '.dat'];
    fid = fopen(line_file,'r');
    mattamp = fscanf(fid,'%f');

    disp('file loaded.')
    fclose(fid);
    inf = mattamp(1:6);
    Lx = inf(1);
    Ly = inf(2);
    nx = inf(3);
    ny = inf(4);
    ntipo=inf(5);
    beta =inf(6);
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
    permmap = zeros(ny,nx);
%
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
%
    if(II==ini)
        vet1=zeros(nx*ny*N,1);
    end
    permmap=Mk*exp(Sk*permmap);
    vd=reshape(permmap,nx*ny,1);
    vet1(num*nx*ny+1:(num+1)*nx*ny,1)=vd;
    num=num+1;
%    media=mean(vd)
%    variancia=var(vd)
%    std=sqrt(variancia)
    clear vd permmap
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
media=mean(vet1)
variancia=var(vet1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('----------------')
disp(' LOADING FIELD ');
%
base_name = [file_base namet2];
%
num=0;
fim=ini+N-1;
for II=ini:1:fim
    line_file = [base_name num2str(II,5) '.dat'];
    fid = fopen(line_file,'r');
    mattamp = fscanf(fid,'%f');

    disp('file loaded.')
    fclose(fid);
    inf = mattamp(1:6);
    Lx = inf(1);
    Ly = inf(2);
    nx = inf(3);
    ny = inf(4);
    ntipo=inf(5);
    beta =inf(6);
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
    permmap = zeros(ny,nx);
%
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
%
    if(II==ini)
        vet2=zeros(nx*ny*N,1);
    end
    permmap=Mp*exp(Sp*permmap);
    vd=reshape(permmap,nx*ny,1);
    vet2(num*nx*ny+1:(num+1)*nx*ny,1)=vd;
    num=num+1;
%    media=mean(vd)
%    variancia=var(vd)
%    std=sqrt(variancia)
    clear vd permmap
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
media=mean(vet2)
variancia=var(vet2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plot(vet1,vet2,'o')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
%close all
