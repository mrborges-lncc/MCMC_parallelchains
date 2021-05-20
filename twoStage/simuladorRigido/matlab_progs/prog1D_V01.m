clear;
close all
Lx = 1.00; % Tamanho do dominio
nx = 500;  % Malha do dominio
dx = Lx/double(nx); % Delta x
T = 0.2687; % Tempo total
phi = 1.0; % porosidade
v = 1.0; % velocidade
srw = 0.0;
sro = 0.0;
s0 = srw;
muw = 1.0;
muo = 20.0;
tfluxo = 1;
%% CONDICAO CFL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
v  = v/phi;
vmax = velocmax(srw,sro,muw,muo,v,tfluxo);
courant = 1;
dt = courant*(dx/vmax);
nt = int16(T/dt)+1;
dt = T/double(nt);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sn = zeros(nx,2)+s0;
s  = zeros(nx,1);
x  = [dx/2.0:dx:Lx-dx/2.0];
xa = [0:Lx/15000:Lx];
sa = zeros(size(xa,2),1)+srw;
%% DEFINICAO DO METODO %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
flag = 1; % Upwind            ---> nn = 2
%flag = 2; % Lax-Friedrichs    ---> nn = 2
%flag = 3; %Leapfrog           ---> nn = 2
%flag = 4; %Lax-Wendroff       ---> nn = 2
%flag = 5; %Beam-Warming       ---> nn = 3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% condicao inicial para sn %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nn = 2;
sn(1:nn-1,1) = 1.0-sro;
if(tfluxo==3)
    saltoi = 0.1;
    saltoj = 0.3;
    for j=1:nx
        if(x(j)>=saltoi&&x(j)<=saltoj)
            sn(j,1)=1-sro;
        else
            sn(j,1)=srw;
        end
    end
end
sn(:,2)=sn(:,1);
s=sn(:,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for k=1:nt-round((x(nn-1)/vmax)/dt)
    for i=nn:nx-2
        if(flag==1)
            %%% Upwind
            s(i,1) = sn(i,1) - (dt/(dx))*(funcF(sn(i,1),v,srw,sro,muw,muo,tfluxo)-funcF(sn(i-1,1),v,srw,sro,muw,muo,tfluxo));
        end
        if(flag==2)
            %%% Lax-Friedrichs
            s(i,1) = 0.5*(sn(i-1,1)+s(i+1,1)) - (0.5*dt/(dx))*(funcF(sn(i+1,1),v,srw,sro,muw,muo,tfluxo)-funcF(sn(i-1,1),v,srw,sro,muw,muo,tfluxo));
        end
        if(flag==3)
            %%%Leapfrog
            s(i,1) = sn(i,2) - (dt)/(dx)*(funcF(sn(i+1,1),v,srw,sro,muw,muo,tfluxo)-funcF(sn(i-1,1),v,srw,sro,muw,muo,tfluxo));     
        end
        if(flag==4)
            %%%Lax-Wendroff
            s(i,1) = sn(i,1) - dt/((2*dx))*(funcF(sn(i+1,1),v,srw,sro,muw,muo,tfluxo)-funcF(sn(i-1,1),v,srw,sro,muw,muo,tfluxo))+((dt^2)*(dfuncF(sn(i,1),srw,sro,muw,muo,v,tfluxo)^2)/(2*(dx^2)))*(funcF(sn(i+1,1),v,srw,sro,muw,muo,tfluxo)-2*(funcF(sn(i,1),v,srw,sro,muw,muo,tfluxo))+funcF(sn(i-1,1),v,srw,sro,muw,muo,tfluxo));     
        end
        if(flag==5)
            %%%Beam-Warming
            s(i,1) = sn(i,1) - (dt/((2*dx)))*(3*(funcF(sn(i,1),v,srw,sro,muw,muo,tfluxo))-(4*(funcF(sn(i-1,1),v,srw,sro,muw,muo,tfluxo)))+funcF(sn(i-2,1),v,srw,sro,muw,muo,tfluxo))+((dt^2)*(dfuncF(sn(i,1),srw,sro,muw,muo,v,tfluxo)^2)/(2*(dx^2)))*(funcF(sn(i,1),v,srw,sro,muw,muo,tfluxo)-(2*(funcF(sn(i-1,1),v,srw,sro,muw,muo,tfluxo)))+funcF(sn(i-2,1),v,srw,sro,muw,muo,tfluxo));   
        end
    end
    sn(:,2)=sn(:,1);
    sn(:,1)=s;
end
    %% solucao analitica %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(tfluxo == 1)
    salto = dfuncF(sn(i,1),srw,sro,muw,muo,v,tfluxo)*T;
    for j=1:size(xa,2)
        if(xa(j)<=salto)
            sa(j,1)=1.0-sro;
        end
    end
end
if(tfluxo == 2)
    sol = load('solanalitica.dat');
    sol = sol';
    sa = sol(:,2);
    xa = sol(:,1);
end
 figc(xa, sa, x, s)
%% Consevação da massa
intanalitica = trapz(xa,sa);
intnumerica  = trapz(x,s);
%% Norma %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[xa, idx] = sort(xa);
y=sa;
for j=1:length(idx)
    sa(j) = y(idx(j));
end
y=interp1(xa,sa,x);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('\n#################################################\n');
fprintf('Integral da solução analítica..........: %2.5f\n',abs(intanalitica));
fprintf('Integral da aproximação numérica.......: %2.5f\n',intnumerica);
fprintf('Norma L^2 entre as curvas..............: %2.5f\n',norm((s-y'),2)/norm(y'));
fprintf('#################################################\n\n');
clear;