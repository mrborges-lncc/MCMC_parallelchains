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
base_aux = '../figuras/';
N=100;
NM=1000;
%[FILENAME, PATHNAME] =uigetfile('../simul_comp/exp/fields/*.dat', 'LOAD DATA');
% [FILENAME, PATHNAME] =uigetfile(...
%     '/Users/mrborges/ESTUDO_GEO/projMCMC/trunk/forecast/3Dfields/fields/*.dat', 'LOAD DATA');
% [FILENAME, PATHNAME] =uigetfile(...
%      '/Users/mrborges/projMCMC/trunk/forecast/3Dfields/fields/MC_0*.dat', 'LOAD DATA');
% % [FILENAME, PATHNAME] =uigetfile('../twoStage/select_fields/*.dat', 'LOAD DATA');
% line_file=sprintf('%s%s', PATHNAME,FILENAME);
line_file='/Users/mrborges/projMCMC/trunk/forecast/3Dfields/fields/MC_0_0.dat'
%line_file='../simuladorBTMM/exp01/fields/ref_0.dat'
%
output ='/Users/mrborges/UW_MCMC/fieldsMC/permMC_';
%
if(N==1)
    n=1;
else
    n=0;
end
for i=size(line_file,2):-1:1
   a=line_file(i);
   if(a=='.')
       in=line_file(i-1);
       for j=i-2:-1:1
           a=line_file(j);
           if(a=='_')
               if(n==0)
                   n=1
               else
                   str_k=j;
                   i=0;
                   break
               end
           else
               in=[a in]
           end
       end
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
base_name = file_base(str_k:end);
%
for files=303:NM
    ini=str2num(in);
    fim=ini+N-1;
    for II=ini:1:fim
        nslice = num2str(II,5);
        nfile  = num2str(files-1,5);
        if(N==1)
            line_file = [file_base nfile '.dat']
        else
            line_file = [file_base nfile '_' nslice '.dat']
        end
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
        Xi=zeros(nx*ny,1);
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
        permmap=[];
        k=0;
        m=0;
        for j=ny:-1:1
            k=k+1;
            if(mattamp(k)~=ny-j)
                disp('erro1')
                break
            end
            for i=1:nx
                k=k+1;
                m=m+1;
                permmap(j,i)=mattamp(k);
                Xi(m)=mattamp(k);
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
    %
        media=mean(mean(permmap))
        vd=reshape(permmap,nx*ny,1);
        variancia=var(vd)
        std=sqrt(variancia)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if(N==1)
            bb2 = [output nfile '.dat']
        else
            bb2 = [output nfile '_' nslice '.dat']
        end
        outfile2 = fopen(bb2, 'wt');
        fprintf(outfile2,'# campos KL\n');
        fprintf(outfile2,'%d x %d x 1\n',nx,ny);
        for i=1:nx*ny
            fprintf(outfile2,'%f\n',Xi(i));
        end      
       fclose(outfile2);
    end
end
!cd /Users/mrborges/UW_MCMC/fieldsMC/; pwd; sh enviar.sh
clear
%close all
