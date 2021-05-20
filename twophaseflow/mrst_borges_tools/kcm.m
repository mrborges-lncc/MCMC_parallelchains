function [A]=kcm(c,phi)
A=c*(phi.^3)./((1-phi));
