function [dw]=dlbw(s,srw,sro,muw,muo)
    lw=2*(s-srw)/(muw*(1-srw)^2);
    lo=-2*(1-sro-s)/(muo*(1-sro)^2);
    lt=lbw(s,srw,muw)+lbo(s,sro,muo);
    dw=lw/(lt)-krw(s,srw)*(lw+lo)/(lt*lt);
