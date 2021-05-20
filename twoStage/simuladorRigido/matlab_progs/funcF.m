function [out] = funcF(s,vel,srw,sro,muw,muo,flag)
 if(flag==1)
     out = s*vel;
 end
 if(flag==2)
     out=lbw(s,srw,muw)/(lbw(s,srw,muw)+lbo(s,sro,muo));
 end
 if(flag==3)
     out = 0.5*s*s;
 end