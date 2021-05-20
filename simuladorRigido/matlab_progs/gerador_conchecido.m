clear
Lx=1.0;
Ly=1.0;
Nx=100;
Ny=100;    
nome = ['barr100x100_0']

%% definicoes da malha %%%%%%%%%%%%%%%%%%%%%%%%%%
num_elem=Nx*Ny;
dx=Lx/Nx;
dy=Ly/Ny;
x=[0+dx*0.5:dx:Lx-dx*0.5];
y=[0+dy*0.5:dy:Ly-dy*0.5];
coord=zeros(2,num_elem);
%% vetor de coordenadas %%%%%%%%%%%%%%%%%%%%%%%%%
k = 0;
for i = 1:Ny
  for j = 1:Nx
    k = k + 1;
    coord(1,k) = (j-1)*dx+dx*0.5;
    coord(2,k) = (i-1)*dy+dy*0.5;
  end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
K=[];
k=0;
n=0.0;
for j=Ny:-1:1
    for i=1:Nx
      n=n+1;
        K(j,i)=log(1.0);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
diam = .02;
px = 0.42;
py = 0.85;
xmax = px+diam*4;
xmin = px-diam*4;
ymax = py+diam;
ymin = py-diam;
%
k=0;
for j=Ny:-1:1
    for i=1:Nx
        k=k+1;
        if(coord(2,k)<=ymax && coord(2,k)>=ymin)
            if(coord(1,k)<=xmax && coord(1,k)>=xmin)
                K(j,i) = log(0.01);
            end
        end
    end
end
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
diam = .02;
px = 0.58;
py = 0.7;
xmax = px+diam*4;
xmin = px-diam*4;
ymax = py+diam;
ymin = py-diam;
%
k=0;
for j=Ny:-1:1
    for i=1:Nx
        k=k+1;
        if(coord(2,k)<=ymax && coord(2,k)>=ymin)
            if(coord(1,k)<=xmax && coord(1,k)>=xmin)
                K(j,i) = log(0.01);
            end
        end
    end
end
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% diam = .2;
% px = 1.5;
% py = 0.75;
% xmax = px+diam;
% xmin = px-diam;
% ymax = py+diam;
% ymin = py-diam;
% %
% k=0;
% for j=Ny:-1:1
%     for i=1:Nx
%         k=k+1;
%         if(coord(2,k)<=ymax && coord(2,k)>=ymin)
%             if(coord(1,k)<=xmax && coord(1,k)>=xmin)
%                 K(j,i) = log(0.5);
%             end
%         end
%     end
% end
% %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Lx=4.0;
% Ly=3.0;
% Nx=4;Ny=3;
% dx=Lx/Nx;
% dy=Ly/Ny;
% x=[0+dx*0.5:dx:Lx-dx*0.5];
% y=[0+dy*0.5:dy:Ly-dy*0.5];
% [X,Y]=meshgrid(x,y)
% K=[1 3 1 1;2 1 3 2; 3 1 2 3]
%K=[1 1 1 1 1;1 1 1 0.1 1;1 1 1 1 1;1 1 1 1 1;1 1 1 1 1]
%     bb = ['/home/mrborges/mixed_hybrid/raviart_maicon/ppg-rt/field/K.dat'];
%     outfile = fopen(bb, 'wt');
%     fprintf(outfile, '   %d %d\n',1,Nx*Ny);
% 
%     for j=Ny:-1:1
%         for i=1:Nx
%             fprintf(outfile,'   %1.7e\t',K(j,i));
%         end
%     end
%     fclose(outfile);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    bb2 = ['../field/' nome '.dat'];
    outfile2 = fopen(bb2, 'wt');
    fprintf(outfile2,'%f\n',Lx);
    fprintf(outfile2,'%f\n',Ly);
    fprintf(outfile2,'%d\n',Nx);
    fprintf(outfile2,'%d\n',Ny);
    fprintf(outfile2,'%d\n',2);
    fprintf(outfile2,'%f\n',0.5);
    fprintf(outfile2,'%d\n',0);
    fprintf(outfile2,'%d\n',0);

%     for j=1:Ny
%       fprintf(outfile2,'%d\n',j-1);
%       for i=1:Nx
%         fprintf(outfile2,'%f ',K(Ny-j+1,i));
%       end
%       fprintf(outfile2,'\n192837465\n');
%     end
    for j=Ny:-1:1
      fprintf(outfile2,'%d\n',Ny-j);
      for i=1:Nx
        fprintf(outfile2,'%f ',K(j,i));
      end
      fprintf(outfile2,'\n192837465\n');
    end
    fclose(outfile2);
display('SALVO')

%     bb3 = ['/home/mrborges/mixed_hybrid/raviart_maicon/ppg-rt-cc/field/rf_50x50_1.dat'];
%     outfile2 = fopen(bb3, 'wt');
%     fprintf(outfile2,'%f\n',Lx);
%     fprintf(outfile2,'%f\n',Ly);
%     fprintf(outfile2,'%d\n',Nx);
%     fprintf(outfile2,'%d\n',Ny);
%     fprintf(outfile2,'%d\n',2);
%     fprintf(outfile2,'%d\n',0);
%     fprintf(outfile2,'%d\n',0);
%     fprintf(outfile2,'%d\n',0);
% %
%     for j=1:Ny
%       fprintf(outfile2,'%d\n',j-1);
%       for i=1:Nx
%         fprintf(outfile2,'%f ',K(j,i));
%         K1(j,i)=K(Ny-j+1,i);
%       end
%       fprintf(outfile2,'\n192837465\n');
%     end
%     fclose(outfile2);
%pcolor(K)
    xelem=[];
    yelem=[];
    elem=[];
    dxx=dx*0.5;
    dyy=dy*0.5;
    figure(1)
    axis equal
    nel=0;
    for j=Ny:-1:1
        for i=1:Nx
            nel=nel+1;
            xelem(1)=coord(1,nel)-dxx;
            xelem(3)=coord(1,nel)+dxx;
            xelem(2)=coord(1,nel)+dxx;
            xelem(4)=coord(1,nel)-dxx;
            yelem(1)=coord(2,nel)-dyy;
            yelem(2)=coord(2,nel)-dyy;
            yelem(3)=coord(2,nel)+dyy;
            yelem(4)=coord(2,nel)+dyy;
            nos(1)=K(j,i);nos(2)=K(j,i);nos(3)=K(j,i);nos(4)=K(j,i);
            patch(xelem,yelem,nos);
            if(nel==1)
                hold on;
            end
        end
    end
            axis equal
    A=Lx;
    B=Ly;
    axis([0 A 0 B]);
    daspect([1 1 1]); 
    base=['../figuras/' nome]
    print('-djpeg90',base)
%    print('-depsc',base)

%clear;