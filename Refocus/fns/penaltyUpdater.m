function [mu_out,mu_update] = penaltyUpdater(mu,primal_res,dual_res,resid_tol,mu_inc,mu_dec)

if primal_res > resid_tol*dual_res
        mu_out = mu*mu_inc;
        mu_update = 1;
    elseif dual_res > resid_tol*primal_res
        mu_out = mu/mu_dec;
        mu_update = -1;
    else
        mu_out = mu;
        mu_update = 0;
end

end