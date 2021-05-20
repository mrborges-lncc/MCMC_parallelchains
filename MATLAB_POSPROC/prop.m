function [y]=prop(m,s)
  y = mvnrnd(m,s,1)';
  %y = c*unifrnd(-10,10,2,1);