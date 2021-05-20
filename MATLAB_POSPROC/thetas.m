function [th] = thetas(x,y,n)
  th = polyfit(x,y,n-1);
end

