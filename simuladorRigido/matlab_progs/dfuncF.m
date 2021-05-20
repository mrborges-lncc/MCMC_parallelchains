function [out] = dfuncF(s,srw,sro,muw,muo,v,flag)
    if(flag == 1)
        out = v;
    end
    if(flag == 2)
        out = dlbw(s,srw,sro,muw,muo);
    end
    if(flag == 3)
        out = s;
    end
