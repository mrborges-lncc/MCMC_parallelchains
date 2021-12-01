function [r]=chooser(j,Nc)
  seq = 1:Nc;
  seq = seq(find(seq ~= j));
  r   = datasample(seq,2,'Replace',false);
  return
end
  