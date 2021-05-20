function [out] = velocmax(srw,sro,muw,muo,v,flag)
    if(flag == 1)
        out = v;
    end
    if(flag == 2 || flag == 3)
        s = [srw:((1-sro)-srw)/1000:1-sro];
        out = -1e30;
        for i=1:length(s)
            if(dfuncF(s(i),srw,sro,muw,muo,v,flag)>=out)
                out = dfuncF(s(i),srw,sro,muw,muo,v,flag);
            end
        end
    end