clear
close all

%% Establishing Rate Parameters (should not really change) and declare molecule names
% t is a struct with fields for each tissue name and values denoting the
% number associated with the tissue.
% m is a struct with fields for each molecule name and values denoting the
%number associated with the molecule.
% p is a struct with fields for each parameter name and values indicating
% the value of the parameter.
% cnames is a 8x1 array with tissue names
% mnames is a 21x1 array with molecules names
% pnmaes is a 71x1 array with parameter names
% initconc is a nested struct with fields for each tissue and values
% containing struct with fields for each moelulce and values denoting the
% initial concentrations.
[c,m,p,cnames,mnames,pnames,initconc] = declareParams_multi_tissue_VEGF();
nm = length(mnames);
nc = length(cnames);

%% Parameter changes
% % test without N1 and N2
% file_appender = 'no_N';
% p = alter_params(file_appender, p);

% test without cross-talk between comparments
file_appender = 'no_perm';
p = alter_params(file_appender, p);

% % Try increased R2 production:
% p.k_prod.R2 = p.k_prod.R2*1e2;
% % Try inreasing k_on for R2-V165
% % p.k_on.R2_V165 = p.k_on.R2_V165*1e10;

%% Initial Conditions 

t_end = 2*60*60;
runVar = 0;

y0=zeros(nm*nc,1);
for i=1:nc
    for j=1:nm
        index = (i-1)*nm + j;
        y0(index)=initconc.(cnames{i}).(mnames{j});
    end
end

%% Run with current V165 secretion value
col_names = cell(nc*nm + 3,1);
col_names{1} = 'q_V165';
col_names{2} = 'kprod_R2';
col_names{3} = 'time';
for(i = 1:nc)
    for(j = 1:nm)
        col_names{(i-1)*nm + j + 3} = strcat(cnames{i}, '.', mnames{j});
    end
end

[t,y] = multi_tissue_main_VEGF(c,p,m,y0,t_end,runVar);
collect_results = array2table([repelem(p.q.V165(c.primary), length(t), 1) repelem(p.k_prod.R2, length(t), 1) t y], ...
    'VariableNames', col_names);

%% Test various values of secretion
% Debug with increased values of R2 and secretion
% secretion = [0.27, 1, 1e1, 1e2];
secretion = [0.027];
prod = p.k_prod.R2.*[1e1 1e2 1e3];
for(q=secretion)
    for(kp=prod)
        p.q.V165 = [q, q, q];
        p.k_prod.R2 = kp;
    % run until t_end:
        [t,y] = multi_tissue_main_VEGF(c,p,m,y0,t_end,runVar);
        results = array2table([repelem(p.q.V165(c.primary), length(t), 1) repelem(p.k_prod.R2, length(t), 1) t y], ...
    'VariableNames', col_names);
        collect_results = vertcat(collect_results, results);
    end
end

writetable(collect_results, strcat("debug_results/simulation_results_debug_R2_", file_appender,".csv"));
% writetable(collect_results, strcat("debug_results/simulation_results_debug_R2.csv"));
