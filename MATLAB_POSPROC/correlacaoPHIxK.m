% permfield.m
% -----------
escala=1;
clear
%prompt2={'Diretorio atual figuras/: '};
%Answers2=inputdlg(prompt2,'Nome base para os arquivos',1);
%base_name = char(Answers2);
%base_aux = '../figuras/';
base_aux = '../figuras/';

N=1;
imprK=10; % se igual a 1 gera o grafico de K

disp('----------------')
disp(' LOADING FIELD ');

%[FILENAME, PATHNAME] = uigetfile('../LU/campos/e*.dat', 'LOAD DATA');
%[FILENAME, PATHNAME] =uigetfile('/home/mrborges/RANDOMFIELDS/trunk/FORTRAN/LABTRANGEO_COND/campos/p*.dat', 'LOAD DATA');
%permeabilidade=sprintf('%s%s', PATHNAME,FILENAME);
permeabilidade='../KL/campos/ec50x50_01_0.dat'
permeabilidade='../../FORTRAN/LABTRANGEO_COND/campos/p512x512_0.dat'
permeabilidade='../twophaseflow/exp/fields/perm_ref_0.dat'
%permeabilidade='../simuladorBTMM/exp01/fields/phiamostra_0.dat'
porosidade='../twophaseflow/exp/fields/poro_ref_0.dat'


M=1.0;
S=1.0;

M2=0.435;
S2=0.265;


n=0;
for i=size(permeabilidade,2):-1:3
    a=permeabilidade(i);
    if(a=='.');
        in=permeabilidade(i-1);
        for j=i-2:-1:1
            a=permeabilidade(j);
            if(a=='_')
                str_k=j;
                i=0;
                break
            else
                in=[a in];
            end
        end
    end
end

file_base=permeabilidade(1:str_k);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n=0;
for i=size(porosidade,2):-1:3
    a=porosidade(i);
    if(a=='.');
        in=porosidade(i-1);
        for j=i-2:-1:1
            a=porosidade(j);
            if(a=='_')
                str_k=j;
                i=0;
                break
            else
                in=[a in];
            end
        end
    end
end
file_base2=porosidade(1:str_k);
ini=str2num(in);
fim=ini+N-1;
k=1;
%
for II=ini:1:fim
    permeabilidade = [file_base num2str(II,5) '.dat']
    fid = fopen(permeabilidade,'r');
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
    if(ini==II)
        pk=zeros(nx*ny*N,1);
        pp=zeros(nx*ny*N,1);
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
            permmap(j,i)=M*exp(S*mattamp(k));
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
    ma=max(max(permmap))
    mi=min(min(permmap))
    med=mean(mean(permmap))

    porosidade = [file_base2 num2str(II,5) '.dat']
    fid2 = fopen(porosidade,'r');
    mattamp = fscanf(fid2,'%f');

    disp('file loaded.')
    fclose(fid2);
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
    phimap=[];
    k=0;
    for j=ny:-1:1
        k=k+1;
        if(mattamp(k)~=ny-j)
            disp('erro1')
            break
        end
        for i=1:nx
            k=k+1;
            phimap(j,i)=M2*exp(S2*mattamp(k));
        end
        k=k+1;
        if(mattamp(k)~=192837465)
            disp('erro2')
            break
        end       
    end
    clear mattamp
    s = size(phimap);
    x = s(:,2);
    y = s(:,1);
    ma=max(max(phimap))
    mi=min(min(phimap))
    med=mean(mean(phimap))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   k=1;
   for i=1:nx
       for j=1:ny
           pk(k)=permmap(j,i);
           pp(k)=phimap(j,i);
           k=k+1;
       end   
   end
end
X=sort(pp,x);
%    plot(log(pp),log(pk),'o')
   plot((pp),(pk),'o')
% 
%    base1=['permg_' base_name]
%    base=[base_aux base1 '_' num2str(II,5)]
%    set(gcf,'PaperPositionMode','auto');
%    print('-depsc','-r1200',base);
%    print('-djpeg99',base);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%clear
