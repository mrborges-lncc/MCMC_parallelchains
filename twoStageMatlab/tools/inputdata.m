function [newexp, expname, prop_method, jump, nStage, num_rockpar, ...
            num_datatype, num_trials, NC, freqj, prt, ...
            file_ref, precision, precision_coarse, data_normal,...
            physical_dim, fine_mesh, coarse_mesh, file_KL, KLM] = ...
            inputdata(test)
    if test
        [newexp, expname, prop_method, jump, nStage, num_rockpar, ...
            num_datatype, num_trials, NC, freqj, prt] = finputbox();
        [file_ref, precision, precision_coarse, data_normal] = ...
            finputbox2(nStage, num_datatype);
        [physical_dim, fine_mesh, coarse_mesh, file_KL, KLM] = ...
            finputbox3(nStage, num_rockpar);
    else
        newexp      = true;
        expname     = teste;
        prop_method = 'CW';
        jump        = 0;
        nStage      = 1;
        num_rockpar = 2;
        num_datatype= 2;
        num_trials  = 5000;
        NC          = 3;
        freqj       = 10;
        prt         = 1;
        file_ref = ['~/MCMC_parallelchains/twophaseflow/exp/pres/pres_ref_0.dat';
            '~/MCMC_parallelchains/twophaseflow/exp/prod/prod_ref_0.dat'];
        precision = [2e-4 2e-5];
        precision_coarse = [2.5e-4 2.5e-5];
        data_normal  = 0;
        physical_dim = [ 510 510 20];
        fine_mesh    = [51 51 5];
        coarse_mesh  = [51 51 5];
        file_KL      = [];
        KLM          = 13005;
    end
end

