function [np jump] = nprint_corr(npr,ndt)
    npr = int16(npr);
    ndt = int16(ndt);
    if npr >= ndt
        np = ndt;
        jump= 1;
    else
        jump = ceil(ndt/npr);
        np   = npr;
    end
end

