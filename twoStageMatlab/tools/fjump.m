function [out] = fjump(jump,method,d,Nc)
    TOL = 1e-07;
    if abs(jump) < TOL
        switch  method
            case 'RW'
                out = (2.38/sqrt(d)).^2;
            case 'CW'
                out = (2.38/sqrt(d));
            case 'DE'
                out = (2.38/sqrt(2*d));
            otherwise
                out = 1.0;
        end
    else
        out = jump;
    end
end

