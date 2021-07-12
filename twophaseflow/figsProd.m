close all
clear all

addpath ./mrst_borges_tools/
addpath ./mrst-2021a/
startup

nome  = 'ref';

wcut = load(['./exp/prod/wcut_' nome '_0.dat']);
oilp = load(['./exp/prod/prod_' nome '_0.dat']);
pres = load(['./exp/pres/presinj_' nome '_0.dat']);
cpres= load(['./exp/pres/pres_' nome '_0.dat']);
satw = load(['./exp/conc/sw_' nome '_0.dat']);
disp = load(['./exp/disp/disp_' nome '_0.dat']);

fig_matrix(satw ,'time ($days$)','$\mathsf{s}_{w}$','p_',234)
fig_matrix(cpres,'time ($days$)','Pressure ($\mathsf{MPa}$)','p_',135)
fig_matrix(wcut ,'time ($days$)','Water cut ($m^3/day$)','W_',335)
fig_matrix(oilp ,'time ($days$)','Oil production ($m^3/day$)','W_',356)
%fig_matrix(disp ,'time ($days$)','$u_{z}$ ($m$)','p_',357)

