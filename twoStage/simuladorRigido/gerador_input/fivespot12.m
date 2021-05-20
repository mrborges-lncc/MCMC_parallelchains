function fivespot12(Lx,Ly,Lz,nx,ny,nz,q,nome)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hx = Lx/double(nx);
hy = Ly/double(ny);
q=q/(1.0*hx);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
numel = nx*ny;
numnp = (nx+1)*(ny+1);
nz = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fid = fopen(nome,'w');
fprintf(fid,'    0\n');
fprintf(fid,'Arquivo de entrada 2D\n');
fprintf(fid,'         1         0         2\n');
fprintf(fid,'%10d%10d%10d%10d\n',numnp,numel,nx,ny);
fprintf(fid,'%10d%10d%10d\n',1,1,0);
fprintf(fid,'%10d%10d%10.4f%10.4f\n',1,4,0.0,0.0);
fprintf(fid,'%30.4f%10.4f\n',Lx,0.0);
fprintf(fid,'%30.4f%10.4f\n',Lx,Ly);
fprintf(fid,'%30.4f%10.4f\n',0.0,Ly);
fprintf(fid,'%10d%10d%10d%10d\n',nx,1,ny,nx+1);
fprintf(fid,'%10d\n',0);
%%% face inferior %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf(fid,'%10d%10d%10d%10d\n',1,nx-1,1,1);
%%% face superior %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf(fid,'%10d%10d%10d%10d\n',...
    (2*nx+1)*(ny)+1,(2*nx+1)*(ny)+nx-1,1,1);
%%% face esquerda %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n1 = nx+1;
n12 = n1+int16(ny/2-1e-5)*(2*nx+1);
n2 = n1+(2*nx+1)*(ny-1);
fprintf(fid,'%10d%10d%10d%10d\n',n1,n12-(2*nx+1),2*nx+1,1);
fprintf(fid,'%10d%10d%10d%10d\n',n12+(2*nx+1),n2,2*nx+1,1);
%%% face direita %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n1 = n1+nx+(2*nx+1);
n2 = n2+nx-(2*nx+1);
fprintf(fid,'%10d%10d%10d%10d\n',n1,n2,2*nx+1,1);
%%% pressao %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%fprintf(fid,'%10d%10d%10d%10d\n',(2*nx+1)*ny,(2*nx+1)*ny,1,0);
%fprintf(fid,'%10d%10d%10d%10d\n',(2*nx+1)*ny+nx,(2*nx+1)*ny+nx,1,0);
%%% face injecao valores %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf(fid,'%10d\n',0);
fprintf(fid,'%10d%10d%10.2e\n%30.2e\n%10d%10d\n',...
    n12,2,q,q,1,0);
%fprintf(fid,'%10d%10d%10.2e\n%30.2e\n%10d%10d\n',...
%    nx+1,2,q*0.25,q*0.25,1,0);
%fprintf(fid,'%10d%10d%10.2e\n%30.2e\n%10d%10d\n',...
%    (2*nx+1)*ny,2,q,q,1,0);
%fprintf(fid,'%10d%10d%10.2e\n%30.2e\n%10d%10d\n',...
%    (2*nx+1)*ny+nx,2,q,q,1,0);
%%% face %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf(fid,'%10d\n',0);
fprintf(fid,'%10d%10d%10d%10d%10d\n',1,1,4,9,0);
fprintf(fid,'%10d%15.2e%10.2e%10.2e\n',1,1.0,1e12,1e12);
fprintf(fid,'%10.1f%10.1f%10.1f\n',0.0,0.0,0.0);
numnp1=(nx+1)*(ny+1);
fprintf(fid,'%10d%10d%10d%10d%10d%10d%10d\n',...
     1,1,1,2,nx+3,nx+2,1);
fprintf(fid,'%10d%10d%10d%10d%10d%10d\n',...
    nx,1,1,ny,nx,nx+1);
fprintf(fid,'%10d\n',0);
fprintf(fid,'%10d%10d%10d%10d%10d%10d%10d\n',...
    1,1,1,nx+2,2*nx+2,nx+1,1);
fprintf(fid,'%10d%10d%10d%10d%10d%10d\n',...
    nx,1,1,ny,nx,2*nx+1);
fprintf(fid,'%10d\n',0);
fprintf(fid,'*end\n');
fclose(fid);
%%% faces do poco de producao %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nomep = [nome(1:end-3) '_prod.dat'];
fid = fopen(nomep,'w');
fprintf(fid,'%10d\n',2);
fprintf(fid,'%10d\n',2);
n = nx*(ny+1)+ny*(nx+1);
fprintf(fid,'%10d\n',n);
fprintf(fid,'%10d\n',n-nx);
%fprintf(fid,'%10d\n',2);
fprintf(fid,'%10d\n',nx);
fprintf(fid,'%10d\n',2*nx+1);
fclose(fid);
