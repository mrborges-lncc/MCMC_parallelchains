clc
Ly=100;
Lx=100;
Ny=50;
dy=Ly/Ny;
y=[dy/2:dy:Ly-dy/2]';
x=[0.;Lx/2+.5];
z=0.5;
worx='MONITOR_WELL_X_'
wory='MONITOR_WELL_Y_'
worz='MONITOR_WELL_Z_'
m=1;
for j=1:length(x)
    for i=1:Ny
        disp([worx num2str(m+6) '= ' num2str(x(j))])
        disp([wory num2str(m+6) '= ' num2str(y(i))])
        disp([worz num2str(m+6) '= ' num2str(z)])
        disp([' '])
        m=m+1;
    end
end
    

% save('monitores.txt','time','-ascii');
% clear