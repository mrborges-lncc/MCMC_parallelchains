function [out] = reverseKlog(perm,beta,rho)
    out = ((log(perm) - log(beta))/rho);
end

