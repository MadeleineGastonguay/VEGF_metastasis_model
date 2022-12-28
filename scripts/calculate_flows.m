function rates = calculate_flows(y,c,p,m,tissue)
%% define number of tissues and molecules
nc = length(fieldnames(c));
mnames = fieldnames(m);
nm = length(mnames);

i = c.(tissue);
tissue_start = (i-1)*nm;

%% Convert k_on, k_c, and secretion to biologically relevant quantities
[k_on, k_c, q] = unit_conversions(p, i);

%% Trafficking rates
%     r_traf_R1         = p.k_prod.R1 - p.k_int.R1      * y(tissue_start + m.R1);
%     r_traf_R2         = p.k_prod.R2 - p.k_int.R2      * y(tissue_start + m.R2);
%     r_traf_N1         = p.k_prod.N1 - p.k_int.N1      * y(tissue_start + m.N1);
%     r_traf_N2         = p.k_prod.N2 - p.k_int.N2      * y(tissue_start + m.N2);

r_traf_R1_V165    =             - p.k_int.R1_V    * y(tissue_start + m.R1_V165);
r_traf_R2_V165    =             - p.k_int.R2_V    * y(tissue_start + m.R2_V165);
r_traf_N1_V165    =             - p.k_int.N1_V    * y(tissue_start + m.N1_V165);
r_traf_N2_V165    =             - p.k_int.N2_V    * y(tissue_start + m.N2_V165);

r_traf_M_V165_R1  =             - p.k_int.M_V_R1  * y(tissue_start + m.Mebm_V165_R1);
r_traf_M_V165_R2  =             - p.k_int.M_V_R2  * y(tissue_start + m.Mebm_V165_R2);

%     r_traf_R1_N1      =             - p.k_int.N1_R1   * y(tissue_start + m.R1_N1);
%     r_traf_R1_N2      =             - p.k_int.N2_R1   * y(tissue_start + m.R1_N2);

r_traf_R2_V165_N1 =             - p.k_int.N1_R2_V * y(tissue_start + m.R2_V165_N1);
r_traf_R2_V165_N2 =             - p.k_int.N2_R2_V * y(tissue_start + m.R2_V165_N2);

%% Secretion
r_V165_sec = q.V165;

%% Transport between tissue and blood compartment
r_perm_V165  = p.k_p(i) * p.gamma * p.S(i) / p.U(i) * (y(m.V165) / p.K_AV(c.blood) - y(tissue_start + m.V165) / p.K_AV(i));
r_lymph_V165 = p.k_l(i) / (p.K_AV(i) * p.U(i)) * y(tissue_start + m.V165);

%% Collect rates
internalization = r_traf_R1_V165 + r_traf_R2_V165;
internalization = internalization + r_traf_N1_V165 + r_traf_N2_V165;
internalization = internalization + r_traf_M_V165_R1 + r_traf_M_V165_R2;
internalization = internalization + r_traf_R2_V165_N1 + r_traf_R2_V165_N2;


rates = table(r_V165_sec, internalization, -r_lymph_V165, -r_perm_V165, ...
    'VariableNames',{'Secretion','Internalization','Lymphatic_drainage', 'Permeability'});

end

