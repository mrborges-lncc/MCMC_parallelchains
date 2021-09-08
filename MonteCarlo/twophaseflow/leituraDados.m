clear all
close all
name = 'ref';
name2= 'amostra';
pres = load(['exp/pres/pres_' name '_0.dat']);
pres2= load(['exp/pres/pres_' name2 '_0.dat']);
plot(pres(:,1),pres(:,2),'-k','LineWidth',2)
hold on
plot(pres2(:,1),pres2(:,2),'ob','LineWidth',1)
pause(2)
close all

prod = load(['./exp/prod/prod_' name '_0.dat']);
prod2= load(['./exp/prod/prod_' name2 '_0.dat']);
plot(prod(:,1),prod(:,2),'-k','LineWidth',2)
hold on
plot(prod(:,1),prod(:,3),'-r','LineWidth',2)
plot(prod(:,1),prod(:,4),'-b','LineWidth',2)
plot(prod(:,1),prod(:,5),'-g','LineWidth',2)
plot(prod2(:,1),prod2(:,2),'ok','LineWidth',1)
plot(prod2(:,1),prod2(:,3),'or','LineWidth',1)
plot(prod2(:,1),prod2(:,4),'ob','LineWidth',1)
plot(prod2(:,1),prod2(:,5),'og','LineWidth',1)
pause(2)
close all

disp =  load(['./exp/disp/disp_' name '_0.dat']);
disp2=  load(['./exp/disp/disp_' name2 '_0.dat']);
plot(disp(:,1),disp(:,2),'-k','LineWidth',2)
hold on
plot(disp(:,1),disp(:,3),'-r','LineWidth',2)
plot(disp(:,1),disp(:,4),'-b','LineWidth',2)
plot(disp(:,1),disp(:,5),'-g','LineWidth',2)
plot(disp2(:,1),disp2(:,2),'ok','LineWidth',1)
plot(disp2(:,1),disp2(:,3),'or','LineWidth',1)
plot(disp2(:,1),disp2(:,4),'ob','LineWidth',1)
plot(disp2(:,1),disp2(:,5),'og','LineWidth',1)
