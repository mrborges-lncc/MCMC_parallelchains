function [ f ] = blackbox(x,a)
  sz = length(a);
  f = zeros(1,length(x));
  for i=1:sz
      f = f+a(i)*x.^(sz-i);
  end
end