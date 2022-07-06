function cc = matern(xv,yv,sigma,eta1,eta2,eta3,nu)
 x = (xv(:,1)-yv(:,1));
 y = (xv(:,2)-yv(:,2));
 z = (xv(:,3)-yv(:,3));
 h = sqrt((x.^2)/(eta1) + (y.^2)/(eta2)+ (z.^2)/(eta3));
 aux = 1.0/gamma(nu)*(2.0^(nu-1));
 aux2= besselk(nu,h);
 if(aux<1.0e-08||aux2>1.e20)
     cc = sigma;
 else
     cc = -sigma*(1.0 - (aux)*((h).^(nu)).*aux2)+sigma;
 end