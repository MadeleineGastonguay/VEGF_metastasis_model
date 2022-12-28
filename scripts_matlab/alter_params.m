function [p] = alter_params(file_appender, p)
% Enumerate parameter options
cases.no_N = 1;
cases.increased_R2 = 2;
cases.no_perm = 3;
cases.new_transport = 4;
cases.no_main = 5;
cases.new_lymph = 6;
cases.only_R2_new_lymph = 7;
cases.no_ligand_new_lymph = 8;

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

    case 4 % Simulations with updated k_l and k_p for primary compartment
        p.k_l(3) = 120/3600; %cm^3/s %Source: Stefanini 2010 (given in mL/hour)
        p.k_p(3) = 4E-7; %Source: Stefanini 2010

    case 5 % Simulations without main compartment
        p.k_p(2) = 0;
        p.k_l(2) = 0;
        p.k_l(3) = 120/3600; %cm^3/s %Source: Stefanini 2010 (given in mL/hour)
        p.k_p(3) = 4E-7; %Source: Stefanini 2010

     case 6 % Simulations with updated k_l and k_p for primary compartment
        p.k_l(3) = 1.79e-5; %cm^3/s %Source: Feilim claculated it
        p.k_p(3) = 4E-7; %Source: Stefanini 2010
        
     case 7 % Simulations with updated k_l and k_p for primary compartment but only R2
        p.k_l(3) = 1.79e-5; %cm^3/s %Source: Feilim claculated it
        p.k_p(3) = 4E-7; %Source: Stefanini 2010
        p.k_prod.R1 = 0;
        p.k_prod.N1 = 0;
        p.k_prod.N2 = 0;
    
    case 8 % Simulations with updated k_l and k_p for primary compartment but no ligand
        p.k_l(3) = 1.79e-5; %cm^3/s %Source: Feilim claculated it
        p.k_p(3) = 4E-7; %Source: Stefanini 2010
        p.q.V165 = zeros(3, 1);

end