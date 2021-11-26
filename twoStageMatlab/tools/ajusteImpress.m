function [ndata,njump] = ajusteImpress(ndata,nstep)
    if ndata > nstep
        ndata = nstep;
    else
        if int64(mod(nstep,ndata)) ~= 0
            ndata = ndata + 1;
            while int64(mod(nstep,ndata)) ~= 0
                ndata = ndata - 1;
            end
        end
    end
    njump = nstep / ndata;
end

