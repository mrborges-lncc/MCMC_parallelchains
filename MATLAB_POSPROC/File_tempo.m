T=500;
NK=50;
time=[0:T/NK:T-T/NK]';
time(1,1)=0.5*T/NK;
save('K_time_change.txt','time','-ascii');
time
size(time)
clear