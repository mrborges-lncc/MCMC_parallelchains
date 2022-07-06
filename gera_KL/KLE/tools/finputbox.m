function [lx,ly,lz,nx,ny,nz,NX,NY,NZ,...
    clx,cly,clz,beta,nu,vari,NR,intp,mm,typec,TIPOINP,filei]= finputbox()
options.Resize='on';
options.WindowStyle='normal';
options.Interpreter='tex';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
prompt = {'L_x:','L_y:','L_z:'};
dlg_title = 'Domain dimentions';
num_lines = [1 40];
def = {'50.0','30.0','30.0'};
answer = inputdlg(prompt,dlg_title,num_lines,def,options);
ans = char(answer);
lx=str2num(ans(1,:));
ly=str2num(ans(2,:));
lz=str2num(ans(3,:));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
prompt = {'n_x:','n_y:','n_z:'};
dlg_title = 'Mesh to compute the covariance matrix';
num_lines = [1 60];
def = {'25','15','15'};
answer = inputdlg(prompt,dlg_title,num_lines,def,options);
ans = char(answer);
nx=str2num(ans(1,:));
ny=str2num(ans(2,:));
nz=str2num(ans(3,:));
nx=int32(nx);
ny=int32(ny);
nz=int32(nz);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
prompt = {'y/n:'};
dlg_title = 'Perform interpolation to a finer mesh?';
num_lines = [1 60];
def = {'n'};
answer = inputdlg(prompt,dlg_title,num_lines,def,options);
ans = char(answer);
NX=nx;
NY=ny;
NZ=nz;
intp=0;
if(ans=='y')
    intp=1;
    prompt = {'N_x:','N_y:','N_z:'};
    dlg_title = 'Field mesh';
    num_lines = [1 40];
    def = {'20','20','20'};
    answer = inputdlg(prompt,dlg_title,num_lines,def,options);
    ans = char(answer);
    NX=str2num(ans(1,:));
    NY=str2num(ans(2,:));
    NZ=str2num(ans(3,:));
%     NX=int32(NX);
%     NY=int32(NY);
%     NZ=int32(NZ);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
prompt = {'Cov. type:'};
list = {'Exponential','Fractal','Squared-Exponential','Matern variogram'};
num_lines = [1 40];
[indx,tf] = listdlg('PromptString',{'Select a covariance type.',''},...
    'SelectionMode','single','ListString',list);
typec=int64(indx);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if typec == 1
    prompt = {'\lambda_x:','\lambda_y:','\lambda_z:'};
    dlg_title = 'Correlation lengths';
    num_lines = [1 40];
    def = {'2.0','2.0','2.0'};
    answer = inputdlg(prompt,dlg_title,num_lines,def,options);
    ans = char(answer);
    clx=str2num(ans(1,:));
    cly=str2num(ans(2,:));
    clz=str2num(ans(3,:));
    fprintf('\n================================\n')
    fprintf('Exponential Covariance:\n l_x = %4.3f\n l_y = %4.3f\n l_z = %4.3f',clx,cly,clz)
    fprintf('\n================================\n')    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
beta = 0.5;
if typec == 2
    prompt = {'\beta:'};
    dlg_title = 'Hurst coeficient';
    num_lines = [1 40];
    def = {'0.50'};
    answer = inputdlg(prompt,dlg_title,num_lines,def,options);
    ans = char(answer);
    beta= str2num(ans(1,:));
    fprintf('\n================================\n')
    fprintf('Fractal Covariance:\n beta = %4.3f',beta)
    fprintf('\n================================\n')    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if typec == 3
    prompt = {'\lambda_x:','\lambda_y:','\lambda_z:'};
    dlg_title = 'Correlation lengths';
    num_lines = [1 40];
    def = {'2.0','2.0','2.0'};
    answer = inputdlg(prompt,dlg_title,num_lines,def,options);
    ans = char(answer);
    clx=str2num(ans(1,:));
    cly=str2num(ans(2,:));
    clz=str2num(ans(3,:));
    fprintf('\n================================\n')
    fprintf('Squared-Exponential Covariance:\n l_x = %4.3f\n l_y = %4.3f\n l_z = %4.3f',clx,cly,clz)
    fprintf('\n================================\n')    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nu = 0.5;
if typec == 4
    prompt = {'\nu:'};
    dlg_title = 'Matern coefficient';
    num_lines = [1 40];
    def = {'0.50'};
    answer = inputdlg(prompt,dlg_title,num_lines,def,options);
    ans = char(answer);
    nu  = str2num(ans(1,:));
    fprintf('\n================================\n')
    fprintf('Matern Variogram:\n beta = %4.3f',nu)
    fprintf('\n================================\n')    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
prompt = {'\sigma^2:'};
dlg_title = 'Variance';
num_lines = [1 40];
def = {'1.0'};
answer = inputdlg(prompt,dlg_title,num_lines,def,options);
ans = char(answer);
vari=str2num(ans(1,:));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
prompt = {'M:'};
dlg_title = 'Number of terms in the serie (if == 0 => maximum)';
num_lines = [1 70];
def = {'0'};
answer = inputdlg(prompt,dlg_title,num_lines,def,options);
ans = char(answer);
mm=str2num(ans(1,:));
%mm=int32(mm);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
prompt = {'NR:'};
dlg_title = 'Number of realizations';
num_lines = [1 40];
def = {'10'};
answer = inputdlg(prompt,dlg_title,num_lines,def,options);
ans = char(answer);
NR=str2num(ans(1,:));
NR=int32(NR);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
prompt = {'y/n:','file:'};
dlg_title = 'Conditioning';
num_lines = [1 40];
def = {'no','./in/input_cond.dat'};
answer = inputdlg(prompt,dlg_title,num_lines,def,options);
ans = char(answer)
if(ans(1,1:3)=='y')
    TIPOINP=1;
else
    TIPOINP=0;
end
filei = ans(2,:)


% tipo_prt  = 1;      % if == 1 imprime campos no formato LNCC cc. formato UW
% paraview_print = 1; % if == 1 imprime arquivo paravizualizacao no paraview
% printa    = 10; % if == 1 salva a matriz T
% printa_cond=10; % if == 1 EXECUTA e salva melhor condicionamento
% TIPOINPUT = 10; % if == 1 entrada dos pontos condicionados (arquivo input_cond.dat)
return

