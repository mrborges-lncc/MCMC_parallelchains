function [schedule1,schedule2,schedule3,schedule4,nstep1,nstep2,...
    nstep3,nstep4,dt] = injectionScheme(bc,W,vinj,Tf,Tinj,Tinjup,Tinjstop,...
    nstep,ndt)
    deltat = Tf/double(nstep);
    TOL    = 1e-07;
    p1     = Tinj;
    p2     = max(0, Tinjup - Tinj);
    if(Tinjstop < TOL)
        p4 = 0.0;
    else
        p4 = Tf - Tinjstop;
    end
    p3     = Tf - (p1 + p2 + p4);
%    p = (p1 + p2 + p3 + p4)/day
    %% First period %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    nstep1 = max(0,ceil(p1/deltat));
    if abs(nstep1) < TOL
        dt1 = [];
        sinal = 0;
        schedule1.step.val     = dt1';
        schedule1.step.control = ones(numel(schedule1.step.val), 1);
        schedule1.control      = struct('W', [], 'bc', []);
    else
        dt1    = p1/double(nstep1) * ones(1,nstep1);
        dt1    = [dt1(1) .* sort(repmat(2.0.^-[1:ndt ndt],1,1)) dt1(2:end)];
        nstep1 = numel(dt1);
        sinal  = 1;
        W(1).val = 0.0;
        schedule1.step.val     = dt1';
        schedule1.step.control = ones(numel(schedule1.step.val), 1);
        schedule1.control      = struct('W', W, 'bc', bc);
    end
    %% Second period %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    nstep2 = max(0,ceil((p2)/deltat));
    if abs(nstep2) < TOL
        dt2 = [];
        schedule2.step.val     = dt2';
        schedule2.step.control = ones(numel(schedule2.step.val), 1);
        schedule2.control      = struct('W', [], 'bc', []);
    else
        dt2 = (p2)/double(nstep2) * ones(1,nstep2);
        if sinal == 0
            dt2   = [dt2(1) .* sort(repmat(2.0.^-[1:ndt ndt],1,1)) dt2(2:end)];
            sinal = 1;
        end
        nstep2 = numel(dt2);
        W(1).val = 0.0;
        schedule2.step.val     = dt2';
        schedule2.step.control = ones(numel(schedule2.step.val), 1);
        schedule2.control      = struct('W', W, 'bc', bc);
    end
    %% Tird period %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    nstep3 = max(0,ceil(p3/deltat));
    if abs(nstep3) < TOL
        dt3 = [];
        schedule3.step.val     = dt3';
        schedule3.step.control = ones(numel(schedule3.step.val), 1);
        schedule3.control      = struct('W', [], 'bc', []);
    else
        dt3 = (p3)/double(nstep3) * ones(1,nstep3);
        if sinal == 0
            dt3   = [dt3(1) .* sort(repmat(2.0.^-[1:ndt ndt],1,1)) dt3(2:end)];
            sinal = 1;
        end
        nstep3 = numel(dt3);
        W(1).val = vinj;
        schedule3.step.val     = dt3';
        schedule3.step.control = ones(numel(schedule3.step.val), 1);
        schedule3.control      = struct('W', W, 'bc', bc);
    end
    %% Forth period %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    nstep4 = max(0,ceil(p4/deltat));
    if abs(nstep4) < TOL
        dt4 = [];
        schedule4.step.val     = dt4';
        schedule4.step.control = ones(numel(schedule4.step.val), 1);
        schedule4.control      = struct('W', [], 'bc', []);
    else
        dt4 = (Tf-Tinjstop)/double(nstep4) * ones(1,nstep4);
        if sinal == 0
            dt4   = [dt4(1) .* sort(repmat(2.0.^-[1:ndt ndt],1,1)) dt4(2:end)];
            sinal = 1;
        end
        nstep4 = numel(dt4);
        W(1).val = 0.0;
        schedule4.step.val     = dt4';
        schedule4.step.control = ones(numel(schedule4.step.val), 1);
        schedule4.control      = struct('W', W, 'bc', bc);
    end
%    nsteps = nstep1+nstep2+nstep3+nstep4
    dt = [dt1 dt2 dt3 dt4];
%    soma = sum(dt)/day
end

