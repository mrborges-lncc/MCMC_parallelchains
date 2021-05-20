
clear;
base_name = 'KL';
M=30;
base_name = [base_name num2str(M,5)]
%file_name = ['../gera_KL/MATLAB/out/aval1_1x1_100x100_l0.2_M10000.dat'];
file_name = ['../gera_KL/MATLAB/out/aval3_1x1_100x100_l0.2_M5000.dat'];
%file_name = ['../gera_KL/MATLAB/out/aval2_1x1_100x100_b0.5_M10000.dat'];
file_name = ['../gera_KL/MATLAB/out/aval3_3.76x3.76x15.4_16x16x64_l0.5x0.5x0.5_M2048.dat'];
file_name = ['../gera_KL/MATLAB/out/aval3_100x100x200_50x50x100_l20x20x600_M6250.dat'];
file_name = ['../gera_KL/MATLAB/out/aval3_100x100x200_20x20x10_l20x20x400_M4000.dat']
dados=load(file_name);
total=sum(dados(1:end));
trun =sum(dados(1:M));
strength = (trun/total);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nome = ['Energia obtida com ' num2str(M,'%8.2d') ' termos na expansao = '...
    num2str(100*strength,'%8.2f')];
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp(nome)
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create figure
figaval(M,dados,strength);
base=['../figuras/aval_' base_name];
set(gcf,'PaperPositionMode','auto');
print('-depsc','-r300',base);
%print('-djpeg90',base)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculo do strength
lim_st = 0.98;
for i=1:size(dados)
    total=sum(dados(1:end));
    trun =sum(dados(1:i));
    strength = (trun/total);
    if(strength>=lim_st)
        M=i;
        break;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nome = ['Numero de termos na expansao KL (M) para obter uma energia de '...
    num2str(strength,'%8.2f') ' = ' num2str(M,'%8.2d')];
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp(nome)
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
clear