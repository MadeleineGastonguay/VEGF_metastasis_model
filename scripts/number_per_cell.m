function n = number_per_cell(conc, tissue, p)
    n = conc*p.N_av*p.ECSA(tissue)/p.ESAV(tissue);
end
