function [ out ] = unifdisc( A,B )
  aux = rand(1,1);
  FA = 0.0;
  FB = 0.0;
  id = 1;
  nc = max(A,B)-min(A,B)+1;
  p  = 1/double(nc);
  for i=1:nc-1
      FA=FA+p;
      FB=FA+p;
      if((aux>=FA) && (aux<FB))
          id = i+1;
      end
  end
  out = id+min(A,B)-1;

end

