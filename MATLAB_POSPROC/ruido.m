function [ output ] = ruido(x,sig)
  s=length(x);
  output = sqrt(sig)*randn(1,s);
end

