clear;
clc
close all
printa = 10;
home = './figuras/newtonSist-';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%x0=[-0.850,-.850].'; %"chute inicial"
%x0=[0.850,.850].'; %"chute inicial"
b= 1.;
TOL=1e-08;
N=200;
et = 2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% x1 = muY; x2 = varY
muK = 4.48169;
varK= 34.51261;
muK = 100;
varK= 100^2;
muY  = log(muK) - 0.5 * log(varK/muK^2 + 1);
varY = log(varK/muK^2 + 1);
x0=[muY,varY].'; %"chute inicial"
dx = 0.5;
dy = 0.5;
xc=[muY,varY].';
%xc=[0.12,-0.270].';
limites = [xc(1)-dx xc(1)+dx xc(2)-dy xc(2)+dy];

syms F f1 f2 f3 x1 x2 x
x = [x1, x2].';
f1(x) = exp(x1 + x2/2) - muK;
f2(x) = exp(2*x1 + x2) * (exp(x2) - 1) - varK;
% f1(x) = log(x1) - 0.5 * log(x2/x1^2 + 1) - muY;
% f2(x) = log(x2/x1^2 + 1) - varY;
F(x) = [f1,f2].';
J(x) = jacobian(F,x);
iJ(x)= inv(J);
Fm = matlabFunction(F);
invJ = matlabFunction(iJ);
i=0;
if printa == 1
    newtonfigura2d(f1,f2,limites,x0)
    base=[home num2str(i)];
    set(gcf,'PaperPositionMode','auto');
    print('-depsc','-r300', base);
end
pause(2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
raizes=[];
erro = [];
F = Fm(x0(1),x0(2));
erro = [erro; i norm(F)];
raizes=[raizes x0];
while(i<N)
    close all
    i=i+1;
    xbar = x0 - invJ(x0(1),x0(2))*F;
    F = Fm(xbar(1),xbar(2));
    erro = [erro; i norm(F)];
    raizes=[raizes xbar];
    if printa == 1
        newtonfigura2d(f1,f2,limites,xbar)
        base=[home num2str(i)];
        set(gcf,'PaperPositionMode','auto');
        print('-depsc','-r300', base);
        pause(2)
    end
    if(norm(F)<TOL)
        break;
    else
        x0 = xbar;
    end
    fprintf('\n Iter: %d-> x=(%6.5f, %6.5f) -> F(x)=(%4.3g, %4.3g)\n',i,xbar(1),xbar(2),Fm(xbar(1),xbar(2)));
end
if(i>=N)
    fprintf('\n%%%%%%%%%%%%%%%%%%%%%%%\n');
    fprintf('\n     Nao convergiu\n')
    fprintf('\n%%%%%%%%%%%%%%%%%%%%%%%\n');
else
    fprintf('\n ########### Convergiu ###########');
    fprintf('\n Iter: %d-> x=(%6.5f, %6.5f) -> F(x)=(%4.3g, %4.3g)\n',i,xbar(1),xbar(2),Fm(xbar(1),xbar(2)));
    fprintf('\n #################################\n');
end