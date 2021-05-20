function [r1,r2]=chooser(j,Nc)
  r1=j;
  r2=j;
  while(r1==j || r2 == j || r1 == r2)
      r1 = unidrnd(Nc);
      r2 = unidrnd(Nc);
  end