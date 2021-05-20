clear;
Lx = 100;
Ly = 100;
nx = 50;
ny = 50;
beta = 0;
ntipo= 3;
name = 'teste_';
N  = 10;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [FILENAME, PATHNAME] =uigetfile('../simul_comp/exp/fields/*.dat',...
%     'LOAD DATA');
[FILENAME, PATHNAME] =uigetfile('../gera_KL/MATLAB/campos/*.dat',...
    'LOAD DATA');
line_file=sprintf('%s%s', PATHNAME,FILENAME);
n=0;
%
for i=size(line_file,2):-1:1
   a=line_file(i);
   if(a=='.')
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
%
file_base=line_file(1:str_k);
home = file_base(1:end-7);

%
iprt=1;
for i=1:N
    line_file = [file_base num2str(i-1,5) '.dat']
    fid = fopen(line_file,'r');
    cab = fscanf(fid,'%c');
    for j=100:-1:3
        st=cab(j-2:j);
        if(st==' x ')
            break
        end
    end
    data= str2num(cab(j+3:end));
    imprime(Lx,Ly,nx,ny,ntipo,beta,data,i,home,name,iprt)
    fclose(fid);
end
clear
%