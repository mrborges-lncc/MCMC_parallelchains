close all
clear all

addpath ./mrst_borges_tools/
addpath ./mrst-2021a/
startup

nome1 = 'ref';
nome2 = 'amostra';

wcut = load(['./exp/prod/wcut_' nome1 '_0.dat']);
oilp = load(['./exp/prod/prod_' nome1 '_0.dat']);
pres = load(['./exp/pres/presinj_' nome1 '_0.dat']);
cpres= load(['./exp/pres/pres_' nome1 '_0.dat']);
satw = load(['./exp/conc/sw_' nome1 '_0.dat']);

wcut2 = load(['./exp000/prod/wcut_' nome2 '_0.dat']);
oilp2 = load(['./exp000/prod/prod_' nome2 '_0.dat']);
pres2 = load(['./exp000/pres/presinj_' nome2 '_0.dat']);
cpres2= load(['./exp000/pres/pres_' nome2 '_0.dat']);
satw2 = load(['./exp000/conc/sw_' nome2 '_0.dat']);

wcut  = [wcut wcut2(:,2:end)];
oilp  = [oilp oilp2(:,2:end)];
pres  = [pres pres2(:,2:end)];
cpres = [cpres cpres2(:,2:end)];
satw  = [satw satw2(:,2:end)];

fig_matrix(satw,'time ($days$)','$\mathsf{s}_{w}$','p_',234)
fig_matrix(cpres,'time ($days$)','Pressure ($\mathsf{MPa}$)','p_',135)
fig_matrix(wcut,'time ($days$)','Water cut ($m^3/day$)','W_',335)
fig_matrix(oilp,'time ($days$)','Oil production ($m^3/day$)','W_',355)

