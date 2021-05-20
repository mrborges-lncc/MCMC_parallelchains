function [ out ] = multinom( nc,p )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
  aux = rand(1,1);
  FA = 0.0;
  FB = 0.0;
  id = 1;
  for i=1:nc-1
      FA=FA+p(i);
      FB=FA+p(i+1);
      if((aux>=FA) && (aux<FB))
          id = i+1;
      end
  end
  out = id;
end

