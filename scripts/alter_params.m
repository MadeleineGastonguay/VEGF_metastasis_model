function [p] = alter_params(file_appender, p)
% Enumerate parameter options
cases.no_N = 1;
cases.increased_R2 = 2;
cases.no_perm = 3;

setting = cases.(file_appender);

switch setting
    case 1 % Simulations with N1 or N2
        p.k_prod.N1 = 0;
        p.k_prod.N2 = 0;

    case 2 % Simulations with increased R2 production
        p.k_prod.R2 = p.k_prod.R2*1e2;

    case 3 % Simulations with no cross-talk between compartments
        p.k_l = zeros(3,1);
        p.k_p = zeros(3,1);

end