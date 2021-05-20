% permfield.m
% -----------
escala=1;
clear
%prompt2={'Diretorio atual figuras/: '};
%Answers2=inputdlg(prompt2,'Nome base para os arquivos',1);
%base_name = char(Answers2);
%base_aux = '../figuras/';
base_aux = '../figuras/';

N=3;
imprK=10; % se igual a 1 gera o grafico de K

disp('----------------')
disp(' LOADING FIELD ');

%[FILENAME, PATHNAME] = uigetfile('../LU/campos/e*.dat', 'LOAD DATA');
%[FILENAME, PATHNAME] =uigetfile('/home/mrborges/RANDOMFIELDS/trunk/FORTRAN/LABTRANGEO_COND/campos/p*.dat', 'LOAD DATA');
%line_file=sprintf('%s%s', PATHNAME,FILENAME);
line_file='../KL/campos/ec50x50_01_0.dat'
line_file='../../FORTRAN/LABTRANGEO_COND/campos/p512x512_0.dat'
line_file='../simuladorBTMM/exp01/fields/ref_0.dat'
%line_file='../simuladorBTMM/exp01/fields/phiamostra_0.dat'
line_file='../twoStage/simuladorBTMM/exp01/fields/amostra_0.dat'
line_file='../twoStage/select_fields/field_KL20_18.dat'
n=0;
for i=size(line_file,2):-1:3
    a=line_file(i);
    if(a=='.');
        in=line_file(i-1);
        for j=i-2:-1:1
            a=line_file(j);
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

file_base=line_file(1:str_k);


M=1.0;
S=1.0;

ini=str2num(in);
fim=ini+N-1;
for II=ini:1:fim
    line_file = [file_base num2str(II,5) '.dat']
    fid = fopen(line_file,'r');
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
            permmap(j,i)=mattamp(k);
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
    figure1 = figure();
    ma=max(max(permmap))
    mi=min(min(permmap))
    med=mean(mean(permmap))

    box('on');
    hold
    
    for j=ny:-1:1
        yi = (ny-j)*dy;
        yf = yi+dy;
        yy = [yi yf; yi yf];
        for i=1:nx
            xf = i*dx;
            xi = xf - dx;
            xx = [xi xi; xf xf];
            p = permmap(j,i)+[0 0; 0 0];
            surf(xx,yy,p);
        end
    end
    view(0,90)
    shading('flat')
    axis('square')
    % shading('interp')
    % shading('faceted')
    %% Create colorbar
    colorbar
% 
    base1=['permg_' base_name]
    base=[base_aux base1 '_' num2str(II,5)]
    set(gcf,'PaperPositionMode','auto');
    print('-depsc','-r1200',base);
    print('-djpeg99',base);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    clear pcolor
end
clear
