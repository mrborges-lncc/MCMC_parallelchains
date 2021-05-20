function [p] = panaliticaRadial(r,rw,q,pw,fluid,perm,h)
    perm = mean(mean(perm));
    p = pw + ((q*fluid.muW(pw))/(2*pi*perm*fluid.rhoWS*h))* ...
        log(r/rw)*fluid.rhoWS;
end

