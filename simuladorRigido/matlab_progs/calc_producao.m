
clear;
disp('----------------')
disp(' LOADING OUTPUT DATA ');

%[FILENAME, PATHNAME] = uigetfile('/home/mrborges/SimuladoresNT/TRACER/sat/*.dat', 'LOAD DATA');
%[FILENAME, PATHNAME] = uigetfile('/home/mrborges/mixed_hybrid/raviart_maicon/ppg-rt-sim-copy/sat/*.dat', 'LOAD DATA');
%[FILENAME, PATHNAME] = uigetfile('/home/mrborges/mixed_hybrid/raviart_maicon/ppg-rt-sim/sat/*.dat', 'LOAD DATA');
[FILENAME, PATHNAME] = uigetfile('../expref/prod/prod_*.dat', 'LOAD DATA');
%[FILENAME, PATHNAME] = uigetfile('/home/mrborges/SIGER/simuladorCO2/out/sa*.res', 'LOAD DATA');
line_file=sprintf('%s%s', PATHNAME,FILENAME);
%line_file='../out/sat200x100_hom.res'
%
% prompt2={'Diretorio atual figuras/: '};
% Answers2=inputdlg(prompt2,'Nome base para os arquivos',1);
% base_name = char(Answers2);
base_aux = '../figuras/prod/';
% base_name=[base_aux base_name];

for i=length(FILENAME):-1:1
    a=FILENAME(i);
    if(a=='.')
        LF=i-1;
    end
end
base_name = [base_aux FILENAME(1:LF)];
line_file=[PATHNAME FILENAME(1:LF-1)];
%line_file= '../prod/prod_fs_hetP_';
name=line_file;
for i=size(line_file,2):-1:1
    a=line_file(i);
    if(a=='_')
        break;
    end
end
name_base=line_file(1:i);
Lx=1.0;
poro=0.010;
Nrand = 1;
ini = 0;
fim = ini+Nrand;
%
MAX=-1;
for i=ini+1:fim
    clear p
    file_name = [name num2str(i-1,5) '.dat'];
    p=load(file_name);
    a=size(p,1);
    if(a>MAX)
        MAX=a;
        I=i;
    end
end
%
prod=zeros(MAX,Nrand);
tempo=zeros(MAX,Nrand);
P=zeros(MAX,Nrand);
for i=ini+1:fim
    clear p
    file_name = [name num2str(i-1,5) '.dat'];
    p=load(file_name);
    m=size(p,1);
    for j=1:m
        tempo(j,i)=p(j,1);
        prod(j,i)=p(j,2);
    end
    t=[0:max(tempo(:,i))/(MAX-1):max(tempo(:,i))]';
    P(:,i)=interp1(tempo(1:m,i),prod(1:m,i),t,'spline');
end

dt=tempo(2)-tempo(1);
clear prod tempo 
y=mean(P,2);
%plot(t,y,'-o');
hold on
A=dt/100;
T=[0:A:max(t)]';
X=zeros(size(T,1),1);
for i=1:size(T,1)
    if(abs(T(i)-Lx*poro)<=A)
        J=i;
        break
    end
end
for i=J:size(T,1)
    X(i)=1.0;
end
%plot(T,X,'-r')
%%% impressao dos mixings %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sinal=0;
for i=size(name_base,2):-1:1
    str=name_base(i);
    if(str=='/')
        sinal=sinal+1;
    end
    if(sinal==2)
        break;
    end
end
%%% impressao dos mixings %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nome_out=[name_base(1:i) 'prod' name_base(i+5:end) 'Nrand' num2str(Nrand,5) '.dat']
fout = fopen(nome_out,'wt');
for j=1:size(t,1)
    fprintf(fout,'%g %g\n',t(j),y(j));
end
fclose(fout);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%            

prodfigure(t, y, T, 1.0-X)%grid on
sinal=0;
for i=size(name_base,2):-1:1
    str=name_base(i);
    if(str=='/')
        sinal=sinal+1;
    end
    if(sinal==2)
        break;
    end
end
name_fig = [name_base(1:i) 'figuras' name_base(i+5:end) num2str(Nrand,5)]

% print('-djpeg99',name_fig);
print('-depsc',name_fig);
MAX=max(y)
MIN=min(y)
hold on
clear