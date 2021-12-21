function [out] = KL(T,theta,numel)
%     tstart =  tic;
    out = zeros(numel,1);
    for el = 1:numel
        out(el) = T(el,:) * theta;
    end
%     telapsed = toc(tstart);
end

