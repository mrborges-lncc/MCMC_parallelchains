clear;
close all;
options.Resize='on';
options.WindowStyle='normal';
options.Interpreter='tex';
TOL=1e-6;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
prompt = {...
    'Nome do experimento:',...
    'Método (RW; DE; DREAM; AM):',...
    'Número de estágios (1 ou 2):',...
    'Número de variáveis de entrada (priors):',...
    'Número máximo de trials:',...
    'Número máximo de variáveis selecionadas:'
    };
%
dlg_title = 'ENTRADA DE DADOS';
num_lines = [1 50];
def = {'blackboxDE','DREAM','1','3','500000','250'};
answer = inputdlg(prompt,dlg_title,num_lines,def,options);
test = int16(size(answer,1));
if(test==0)
    clear all
    %break
end
ans = char(answer);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nome_exp = char(ans(1,:));
metodo   = char(ans(2,:));
nstagio  = int32(str2num(ans(3,:)));
npriors  = int32(str2num(ans(4,:)));
nmaxtrial= int32(str2num(ans(5,:)));
nmaxselec= int32(str2num(ans(6,:)));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
scup = zeros(npriors,1);
if(nstagio>1)
    for i=1:npriors
        up = ['Upscaling método variável ' num2str(i,3) ': '];
        prompt = {...
        up
        };
    dlg_title = 'Upscaling';
    num_lines = [1 50];
    def = {num2str(i-1,3)};
    answer = inputdlg(prompt,dlg_title,num_lines,def,options);
    test = int16(size(answer,1));
    if(test==0)
        clear all
        %break
    end
    ans = char(answer);
    scup(i,1) = int16(str2num(ans));
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
prompt = {...
    'Quantidade de TIPOS de dados observados:',...
    'Dados serão normalizados? (1-sim; 0-não):',...
    'Simulador usado:'
    };
%
dlg_title = 'ENTRADA DE DADOS';
num_lines = [1 50];
def = {'3','1','0','3'};
answer = inputdlg(prompt,dlg_title,num_lines,def,options);
test = int16(size(answer,1));
if(test==0)
    clear all
    %break
end
ans = char(answer);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ndados   = int16(str2num(ans(1,:)));
normaliz = int16(str2num(ans(2,:)));
nsimul   = int16(str2num(ans(3,:)));
set = zeros(ndados,4);
for i=1:ndados
    prompt = {...
        'Número de pontos de observação:',...
        'Tipo de likelihood:',...
        'Número da avaliação inicial:',...
        'Número da avaliação final..:'
        };
    %
    dlg_title = ['ENTRADA DE DADOS - Conjunto ' num2str(i)];
    num_lines = [1 50];
    def = {'1','2','1','21'};
    answer = inputdlg(prompt,dlg_title,num_lines,def,options);
    test = int16(size(answer,1));
    if(test==0)
        clear all
        %break
    end
    ans = char(answer);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    set(i,1) = int16(str2num(ans(1,:)));
    set(i,2) = int16(str2num(ans(2,:)));
    set(i,3) = int16(str2num(ans(3,:)));
    set(i,4) = int16(str2num(ans(4,:)));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ndados   = int16(str2num(ans(1,:)));
normaliz = int16(str2num(ans(2,:)));
nsimul   = int16(str2num(ans(3,:)));
set2     = zeros(npriors,4);
for i=1:npriors
    prompt = {...
        'Gerador do campo aleatório...:',...
        'Método (opção de proposal)...:',...
        'Dimensão estocástica.........:',...
        'Parâmetro do jump (salto)....:'
        };
    %
    dlg_title = ['ENTRADA DE DADOS - Conjunto ' num2str(i)];
    num_lines = [1 50];
    def = {'6','2','4','1.19'};
    answer = inputdlg(prompt,dlg_title,num_lines,def,options);
    test = int16(size(answer,1));
    if(test==0)
        clear all
        %break
    end
    ans = char(answer);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    set2(i,1) = int16(str2num(ans(1,:)));
    set2(i,2) = int16(str2num(ans(2,:)));
    set2(i,3) = int16(str2num(ans(3,:)));
    set2(i,4) = double(str2num(ans(4,:)));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
errosc=zeros(npriors,1);
errosf=zeros(npriors,1);
for i=1:npriors
    prompt = {...
        ['valor:']
        };
    %
    dlg_title = ['Erro para a variável ' num2str(i)];
    num_lines = [1 50];
    def = {'9.0e-03'};
    answer = inputdlg(prompt,dlg_title,num_lines,def,options);
    test = int16(size(answer,1));
    if(test==0)
        clear all
        %break
    end
    ans = char(answer);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    errosf(i,1) =  double(str2num(ans(1,:)));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(nstagio==2)
    for i=1:npriors
        prompt = {...
            ['valor:']
            };
        %
        dlg_title = ['Erro para a variável ' num2str(i)];
        num_lines = [1 50];
        def = {'9.0e-03'};
        answer = inputdlg(prompt,dlg_title,num_lines,def,options);
        test = int16(size(answer,1));
        if(test==0)
            clear all
            %break
        end
        ans = char(answer);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        errosc(i,1) =  double(str2num(ans(1,:)));
    end
else
    errosc=errosf;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Arquivo de entrada %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
file = ['../twoStage/in/entra_' nome_exp '.in']
fid = fopen(file,'w');
fprintf(fid,'%10d\n',nstagio);
fprintf(fid,'%10d\n',npriors);
st = 'amostra'
for i=1:npriors
    fprintf(fid,'%s\n',st);
end
for i=1:npriors
    fprintf(fid,'%10d',scup(i,1));
end
fprintf(fid,'\n')
fprintf(fid,'%10d\n',nmaxtrial);
fprintf(fid,'%10d\n',nmaxselec);
fprintf(fid,'%10d\n',nsimul);
fprintf(fid,'%10d%10d\n',npriors,normaliz);
for i=1:npriors
    fprintf(fid,'%10d%10d%10d%10d\n',set(i,:));
end
for i=1:npriors
    fprintf(fid,'%10d%10d%10d%10d%12.5e\n',set2(i,1),set2(i,2),0,...
        set2(i,3),set2(i,4));
end
for i=1:npriors
    fprintf(fid,'%12.5e\n',errosf(i,1));
end
for i=1:npriors
    fprintf(fid,'%12.5e\n',errosc(i,1));
end
fclose(fid);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% x   =double(xx);
% y   =double(yy);
% ponto_inicial1 =[xx yy];
% Lx  =double(str2num(ans(5,:)));
% Ly  =double(str2num(ans(6,:)));
% nx  =int16(str2num(ans(7,:))+TOL);
% ny  =int16(str2num(ans(8,:))+TOL);
% flux=double(str2num(ans(9,:)));
% flux=flux/Ly;
% coluna =double(str2num(ans(10,:)));
% tinit  =double(str2num(ans(11,:)));
% ttotal =double(str2num(ans(12,:)));
% ncalc  =int16(str2num(ans(13,:))+TOL);
% tinj   =double(str2num(ans(14,:)));
% pref   =double(str2num(ans(15,:)));
% cotaref=double(str2num(ans(16,:)));
% BHP    =double(str2num(ans(17,:)));
% BHC    =double(yy);
% tdescomp=double(str2num(ans(18,:)));
% malhareg=int16(str2num(ans(19,:))+TOL);
% solver  =char(answer(20,:));
% domo_flag=0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
