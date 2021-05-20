function slaby(Lx,Ly,Lz,nx,ny,nz,q,nome)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
q=q/(Lx*Lz);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hx = Lx/double(nx);
hy = Ly/double(ny);
hz = Lz/double(nz);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
numel = nx*ny*nz;
numnp = (nx+1)*(ny+1)*(nz+1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fid = fopen(nome,'w');
fprintf(fid,'    0\n');
fprintf(fid,'Arquivo de entrada 3D\n');
fprintf(fid,'         1         0         3         1\n');
fprintf(fid,'%10d%10d%10d%10d%10d\n',numnp,numel,nx,ny,nz);
fprintf(fid,'%10d%10d%10d\n',1,1,1);
fprintf(fid,'%10d%10d%10.4f%10.4f%10.4f\n',1,8,0.0,0.0,0.0);
fprintf(fid,'%30.4f%10.4f%10.4f\n',Lx,0.0,0.0);
fprintf(fid,'%30.4f%10.4f%10.4f\n',Lx,Ly,0.0);
fprintf(fid,'%30.4f%10.4f%10.4f\n',0.0,Ly,0.0);
fprintf(fid,'%30.4f%10.4f%10.4f\n',0.0,0.0,Lz);
fprintf(fid,'%30.4f%10.4f%10.4f\n',Lx,0.0,Lz);
fprintf(fid,'%30.4f%10.4f%10.4f\n',Lx,Ly,Lz);
fprintf(fid,'%30.4f%10.4f%10.4f\n',0.0,Ly,Lz);
fprintf(fid,'%10d%10d%10d%10d%10d%10d\n',nx,1,ny,nx+1,nz,(nx+1)*(ny+1));
fprintf(fid,'%10d\n',0);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf(fid,'%10d%10d%10d%10d\n',1,nx*ny,1,1);
fprintf(fid,'%10d%10d%10d%10d\n',...
    nx*ny*nz+((nx+1)*ny+(ny+1)*nx)*nz,...
    1*nx*ny*nz+((nx+1)*ny+(ny+1)*nx)*nz+nx*ny,1,1);
%%% face esquerda %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n1 = nx*ny+1;
n2 = nx*ny+nx;
fprintf(fid,'%10d%10d%10d%10d\n',n1,n2,1,1);
for i=1:nz-1
   fprintf(fid,'%10d%10d%10d%10d\n',...
       i*nx*ny+n1+((nx+1)*ny+(ny+1)*nx)*i,...
       i*nx*ny+n2+((nx+1)*ny+(ny+1)*nx)*i,1,1);
end
%%% face %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n1 = n1+(nx+1)*ny+(ny)*nx;
n2 = n2+(nx+1)*ny+(ny)*nx;
% fprintf(fid,'%10d%10d%10d%10d\n',n1,n2,1,1);
% for i=1:nz-1
%    fprintf(fid,'%10d%10d%10d%10d\n',...
%        i*nx*ny+n1+((nx+1)*ny+(ny+1)*nx)*i,...
%        i*nx*ny+n2+((nx+1)*ny+(ny+1)*nx)*i,1,1);
% end
%%% face %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n1 = nx*ny+nx+1;
n2 = (n1)+(2*nx+1)*(ny-1);
for i=1:nz
   fprintf(fid,'%10d%10d%10d%10d\n',...
       n1+((nx+1)*ny+(ny+1)*nx)*(i-1)+(i-1)*nx*ny,...
       n2+((nx+1)*ny+(ny+1)*nx)*(i-1)+(i-1)*nx*ny,2*nx+1,1);
end
%%% face %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n1 = nx*ny+2*nx+1;
n2 = (n1)+(2*nx+1)*(ny-1);
for i=1:nz
   fprintf(fid,'%10d%10d%10d%10d\n',...
       n1+((nx+1)*ny+(ny+1)*nx)*(i-1)+(i-1)*nx*ny,...
       n2+((nx+1)*ny+(ny+1)*nx)*(i-1)+(i-1)*nx*ny,2*nx+1,1);
end
%%% face injecao valores %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n1 = nx*ny+1;
n2 = nx*ny+nx;
fprintf(fid,'%10d\n',0);
fprintf(fid,'%10d%10d%10.2e\n%30.2e\n%10d%10d\n',n1,2,q,q,nx-1,1);
for i=1:nz-1
    fprintf(fid,'%10d%10d%10.2e\n%30.2e\n%10d%10d\n',...
        i*nx*ny+n1+((nx+1)*ny+(ny+1)*nx)*i,2,q,q,nx-1,1);
end
%%% face injecao valores %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf(fid,'%10d\n',0);
fprintf(fid,'%10d%10d%10d%10d%10d\n',1,1,8,8,0);
fprintf(fid,'%10d%15.2e%10.2e%10.2e\n',1,1.0,1e8,1e10);
fprintf(fid,'%10.1f%10.1f%10.1f\n',0.0,0.0,0.0);
numnp1=(nx+1)*(ny+1);
fprintf(fid,'%10d%10d%10d%10d%10d%10d%10d%10d%10d%10d%10d\n',...
     1,1,1,2,nx+3,nx+2,numnp1+1,numnp1+2,numnp1+nx+3,numnp1+nx+2,1);
fprintf(fid,'%10d%10d%10d%10d%10d%10d%10d%10d%10d\n',...
    nx,1,1,ny,nx,nx+1,nz,nx*ny,numnp1);
fprintf(fid,'%10d\n',0);
n1 = nx*ny+1;
n2 = nx*ny+nx;
n1 = n1+(nx+1)*ny+(ny)*nx;
n2 = n2+(nx+1)*ny+(ny)*nx;
fprintf(fid,'%10d%10d%10d%10d%10d%10d%10d%10d%10d\n',...
    1,1,(nx*ny)+1,nx*ny+1+nx+1,nx*ny+1+nx+1+nx,nx*ny+1+nx,1,...
    n2+1,1);
fprintf(fid,'%10d\n',0);
fprintf(fid,'*end\n');
fclose(fid);
%%% faces do poco de producao %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nomep = [nome(1:end-3) '_prod.dat'];
fid = fopen(nomep,'w');
n1 = nx*ny+1;
n2 = nx*ny+nx;
fprintf(fid,'%10d\n',1);
fprintf(fid,'%10d\n',ny*nz);

for j=n1:n2
    fprintf(fid,'%10d\n',j);
end
for i=1:nz-1
    ni = i*nx*ny+n1+((nx+1)*ny+(ny+1)*nx)*i;
    nf = i*nx*ny+n2+((nx+1)*ny+(ny+1)*nx)*i;
    np = 1;
    for j=ni:np:nf
        fprintf(fid,'%10d\n',j);
    end
end
fclose(fid);
