function [ f ] = blackbox(x,a)
  sz = length(a);
  f = zeros(length(x),1);
  for i=sz:-1:1
      f = f+a(i)*x.^(sz-i);
  end
end