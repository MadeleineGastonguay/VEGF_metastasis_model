function [k_on, k_c, q] = unit_conversions(p, tissue)

    %% Convert k_ons and k_cs to (moles/cm^3 tissue) ^-1 s^-1
    knames = fieldnames(p.k_on);
    for j=1:9  % original k_on in M^-1s^-1
        k_on.(knames{j}) = p.k_on.(knames{j}) * 1000 / p.K_AV(tissue);
    end
    for j=10:13 % original k_on in (moles/cm^2)^-1s^-1
        k_on.(knames{j}) = p.k_on.(knames{j}) / p.ESAV(tissue);
    end
    % original k_c in (moles/cm^2)^-1s^-1
    k_c.N1_R1 = p.k_c.N1_R1 / p.ESAV(tissue);
    k_c.N2_R1 = p.k_c.N2_R1 / p.ESAV(tissue);
    
    %% Convert ligad secretion to moles/cm^3
    qnames = fieldnames(p.q);
    for j=1:length(qnames)
        q.(qnames{j}) = p.q.(qnames{j})(tissue) * p.PSAV(tissue) / (1e-5 * p.N_av);
    end
end