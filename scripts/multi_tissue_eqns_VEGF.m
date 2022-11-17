function dydt = multi_tissue_eqns_VEGF(t,y,c,p,m)
%% define number of tissues and molecules
nc = length(fieldnames(c));
mnames = fieldnames(m);
nm = length(mnames);

dydt = zeros(nc*nm,1);


for i=2:nc %iterate through tissues, excluding blood
    tissue_start = (i-1)*nm;
    
    %% Convert k_ons and k_cs to (moles/cm^3 tissue) ^-1 s^-1
    knames = fieldnames(p.k_on);
    for j=1:9  % original k_on in M^-1s^-1
        k_on.(knames{j}) = p.k_on.(knames{j}) * 1000 / p.K_AV(i);
    end
    for j=10:13 % original k_on in (moles/cm^2)^-1s^-1
        k_on.(knames{j}) = p.k_on.(knames{j}) / p.ESAV(i);
    end
    % original k_c in (moles/cm^2)^-1s^-1
    k_c.N1_R1 = p.k_c.N1_R1 / p.ESAV(i);
    k_c.N2_R1 = p.k_c.N2_R1 / p.ESAV(i);
    
    %% Convert ligad secretion to moles/cm^3
    q = p.q;
    qnames = fieldnames(p.q);
    for j=1:length(qnames)
        q.(qnames{j}) = p.q.(qnames{j})(i) * p.PSAV(i) / (1e-5 * p.N_av);
    end

    %% Binding rates
    r_bind_Mecm_V165   = k_on.M_V165     * y(tissue_start + m.Mecm)    * y(tissue_start + m.V165) - p.k_off.M_V165     * y(tissue_start + m.Mecm_V165);
    r_bind_Mebm_V165   = k_on.M_V165     * y(tissue_start + m.Mebm)    * y(tissue_start + m.V165) - p.k_off.M_V165     * y(tissue_start + m.Mebm_V165);
    r_bind_Mpbm_V165   = k_on.M_V165     * y(tissue_start + m.Mpbm)    * y(tissue_start + m.V165) - p.k_off.M_V165     * y(tissue_start + m.Mpbm_V165);  

    r_bind_R1_V165     = k_on.R1_V165    * y(tissue_start + m.R1)      * y(tissue_start + m.V165) - p.k_off.R1_V165    * y(tissue_start + m.R1_V165);
    r_bind_R2_V165     = k_on.R2_V165    * y(tissue_start + m.R2)      * y(tissue_start + m.V165) - p.k_off.R2_V165    * y(tissue_start + m.R2_V165);
    r_bind_N1_V165     = k_on.N1_V165    * y(tissue_start + m.N1)      * y(tissue_start + m.V165) - p.k_off.N1_V165    * y(tissue_start + m.N1_V165);
    r_bind_N2_V165     = k_on.N2_V165    * y(tissue_start + m.N2)      * y(tissue_start + m.V165) - p.k_off.N2_V165    * y(tissue_start + m.N2_V165);

    r_bind_N1_V165_R2  = k_on.N1_V165_R2 * y(tissue_start + m.N1_V165) * y(tissue_start + m.R2)   - p.k_off.N1_V165_R2 * y(tissue_start + m.R2_V165_N1);
    r_bind_R2_V165_N1  = k_on.R2_V165_N1 * y(tissue_start + m.R2_V165) * y(tissue_start + m.N1)   - p.k_off.R2_V165_N1 * y(tissue_start + m.R2_V165_N1);

    r_bind_N2_V165_R2  = k_on.N2_V165_R2 * y(tissue_start + m.N2_V165) * y(tissue_start + m.R2)   - p.k_off.N2_V165_R2 * y(tissue_start + m.R2_V165_N2);
    r_bind_R2_V165_N2  = k_on.R2_V165_N2 * y(tissue_start + m.R2_V165) * y(tissue_start + m.N2)   - p.k_off.R2_V165_N2 * y(tissue_start + m.R2_V165_N2);
    
    r_bind_R1_N1       = k_c.N1_R1       * y(tissue_start + m.N1)      * y(tissue_start + m.R1)   - p.k_off.N1_R1      * y(tissue_start + m.R1_N1);
    r_bind_R1_N2       = k_c.N2_R1       * y(tissue_start + m.N2)      * y(tissue_start + m.R1)   - p.k_off.N2_R1      * y(tissue_start + m.R1_N2);

    % Formation of MVR complexes:
    % Matrix-bound V165 binding to R1 and R2:
    r_bind_Mebm_V165_R1 = k_on.M_V165_R1 * p.f(i) * y(tissue_start + m.Mebm_V165) * y(tissue_start + m.R1)   - p.k_off.M_V165_R1 * y(tissue_start + m.Mebm_V165_R1);
    r_bind_Mebm_V165_R2 = k_on.M_V165_R2 * p.f(i) * y(tissue_start + m.Mebm_V165) * y(tissue_start + m.R2)   - p.k_off.M_V165_R2 * y(tissue_start + m.Mebm_V165_R2);
    % Matrix binding to R-V165 complex:
    r_bind_R1_V165_Mebm = k_on.R1_V165_M * p.f(i) * y(tissue_start + m.R1_V165)   * y(tissue_start + m.Mebm) - p.k_off.R1_V165_M * y(tissue_start + m.Mebm_V165_R1);  
    r_bind_R2_V165_Mebm = k_on.R2_V165_M * p.f(i) * y(tissue_start + m.R2_V165)   * y(tissue_start + m.Mebm) - p.k_off.R2_V165_M * y(tissue_start + m.Mebm_V165_R2);  
    
    %% Trafficking rates
    r_traf_R1         = p.k_prod.R1 - p.k_int.R1      * y(tissue_start + m.R1);
    r_traf_R2         = p.k_prod.R2 - p.k_int.R2      * y(tissue_start + m.R2);
    r_traf_N1         = p.k_prod.N1 - p.k_int.N1      * y(tissue_start + m.N1);
    r_traf_N2         = p.k_prod.N2 - p.k_int.N2      * y(tissue_start + m.N2);

    r_traf_R1_V165    =             - p.k_int.R1_V    * y(tissue_start + m.R1_V165);
    r_traf_R2_V165    =             - p.k_int.R2_V    * y(tissue_start + m.R2_V165);
    r_traf_N1_V165    =             - p.k_int.N1_V    * y(tissue_start + m.N1_V165);
    r_traf_N2_V165    =             - p.k_int.N2_V    * y(tissue_start + m.N2_V165);

    r_traf_M_V165_R1  =             - p.k_int.M_V_R1  * y(tissue_start + m.Mebm_V165_R1);
    r_traf_M_V165_R2  =             - p.k_int.M_V_R2  * y(tissue_start + m.Mebm_V165_R2);

    r_traf_R1_N1      =             - p.k_int.N1_R1   * y(tissue_start + m.R1_N1);
    r_traf_R1_N2      =             - p.k_int.N2_R1   * y(tissue_start + m.R1_N2);

    r_traf_R2_V165_N1 =             - p.k_int.N1_R2_V * y(tissue_start + m.R2_V165_N1);
    r_traf_R2_V165_N2 =             - p.k_int.N2_R2_V * y(tissue_start + m.R2_V165_N2);

    %% Secretion
    r_V165_sec = q.V165;

    %% Transport between tissue and blood compartment
    r_perm_V165  = p.k_p(i) * p.gamma * p.S(i) / p.U(i) * (y(m.V165) / p.K_AV(c.blood) - y(tissue_start + m.V165) / p.K_AV(i));
    r_lymph_V165 = p.k_l(i) / (p.K_AV(i) * p.U(i)) * y(tissue_start + m.V165);  
        
    %% Equations
    
    dydt(tissue_start + m.V165) = r_V165_sec - r_bind_Mecm_V165 - r_bind_Mebm_V165 - r_bind_Mpbm_V165;
    dydt(tissue_start + m.V165) = dydt(tissue_start + m.V165) - r_bind_R1_V165 - r_bind_R2_V165;
    dydt(tissue_start + m.V165) = dydt(tissue_start + m.V165) - r_bind_N1_V165 - r_bind_N2_V165;
    dydt(tissue_start + m.V165) = dydt(tissue_start + m.V165) + r_perm_V165 - r_lymph_V165;
    
    dydt(tissue_start + m.Mecm) = - r_bind_Mecm_V165;
    dydt(tissue_start + m.Mebm) = - r_bind_Mebm_V165 - r_bind_R1_V165_Mebm - r_bind_R2_V165_Mebm;
    dydt(tissue_start + m.Mpbm) = - r_bind_Mpbm_V165;

    dydt(tissue_start + m.R1) = r_traf_R1 - r_bind_R1_V165 - r_bind_R1_N1 - r_bind_R1_N2 - r_bind_Mebm_V165_R1;
    dydt(tissue_start + m.R2) = r_traf_R2 - r_bind_R2_V165 - r_bind_N1_V165_R2 - r_bind_N2_V165_R2 - r_bind_Mebm_V165_R2;
    dydt(tissue_start + m.N1) = r_traf_N1 - r_bind_N1_V165 - r_bind_R1_N1 - r_bind_R2_V165_N1;
    dydt(tissue_start + m.N2) = r_traf_N2 - r_bind_N2_V165 - r_bind_R1_N2 - r_bind_R2_V165_N2;

    dydt(tissue_start + m.R1_N1) = r_traf_R1_N1 + r_bind_R1_N1;
    dydt(tissue_start + m.R1_N2) = r_traf_R1_N2 + r_bind_R1_N2;

    dydt(tissue_start + m.Mecm_V165) = r_bind_Mecm_V165;
    dydt(tissue_start + m.Mebm_V165) = r_bind_Mebm_V165 - r_bind_Mebm_V165_R1 - r_bind_Mebm_V165_R2;
    dydt(tissue_start + m.Mpbm_V165) = r_bind_Mpbm_V165;

    dydt(tissue_start + m.R1_V165) = r_traf_R1_V165 + r_bind_R1_V165 - r_bind_R1_V165_Mebm;
    dydt(tissue_start + m.R2_V165) = r_traf_R2_V165 + r_bind_R2_V165 - r_bind_R2_V165_N1 - r_bind_R2_V165_N2 - r_bind_R2_V165_Mebm;
    dydt(tissue_start + m.N1_V165) = r_traf_N1_V165 + r_bind_N1_V165 - r_bind_N1_V165_R2;
    dydt(tissue_start + m.N2_V165) = r_traf_N2_V165 + r_bind_N2_V165 - r_bind_N2_V165_R2;

    dydt(tissue_start + m.R2_V165_N1) = r_traf_R2_V165_N1 + r_bind_R2_V165_N1 + r_bind_N1_V165_R2;
    dydt(tissue_start + m.R2_V165_N2) = r_traf_R2_V165_N2 + r_bind_R2_V165_N2 + r_bind_N2_V165_R2;

    dydt(tissue_start + m.Mebm_V165_R1) = r_traf_M_V165_R1 + r_bind_Mebm_V165_R1 + r_bind_R1_V165_Mebm;
    dydt(tissue_start + m.Mebm_V165_R2) = r_traf_M_V165_R2 + r_bind_Mebm_V165_R2 + r_bind_R2_V165_Mebm;
 end

 %% Equations in blood
 
% Collect lymphatic drainage and permeability from each tissue for blood
% equations
r_perm_V165  = 0;
r_lymph_V165 = 0;

for i = 2:nc
    tissue_start = (i-1)*nm;
    r_perm_V165  = r_perm_V165  + p.k_p(i) * p.gamma * p.S(i) / p.U(c.blood) * (y(tissue_start + m.V165) / p.K_AV(i)  - y(m.V165) / p.K_AV(c.blood));
    r_lymph_V165 = r_lymph_V165 + p.k_l(i)                    / p.U(c.blood) *  y(tissue_start + m.V165) / p.K_AV(i);
end

 dydt(m.V165) = - p.k_cl * y(m.V165) + r_perm_V165 + r_lymph_V165;


end

