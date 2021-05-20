function [a]=nint(val)
    if(abs(val)<=1e-3)
        a=0.0;
    else
        i=int32(val);
        aux=abs(val-double(i));
        if(aux>=0.5)
            a=double(i);
        else
            if(val>0)
                a=double(i+1);
            else
                a=double(i-1);
            end
        end
    end
    
    