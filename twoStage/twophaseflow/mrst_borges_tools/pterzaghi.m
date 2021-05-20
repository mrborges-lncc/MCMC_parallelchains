function [p] = pterzaghi(z,t,p0,L,c)
    N = numel(z);
    p = zeros(N,1);
    aux = (4/pi) * p0;
    for k=1:N
        zz = z(k);
        sum= 0.0;
        for m=0:20000
            sum = sum + ((1) / (2*m + 1)) *...
                exp(-(((2*m + 1)^2) * pi * pi * c * t) / ...
                (4 * L^2)) * sin(((2*m + 1)*pi*zz) / (2*L));
        end
        p(k) = aux*sum;
    end
end

