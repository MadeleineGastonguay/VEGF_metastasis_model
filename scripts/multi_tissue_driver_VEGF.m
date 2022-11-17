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

[t,y] = multi_tissue_main_VEGF(c,p,m,y0,t_end,runVar);

%% Write out results to a csv file
col_names = cell(nc*nm + 2,1);
col_names{1} = 'q_V165';
col_names{2} = 'time';
for(i = 1:nc)
    for(j = 1:nm)
        col_names{(i-1)*nm + j + 2} = strcat(cnames{i}, '.', mnames{j});
    end
end

time = repmat(t, nc, 1);
compartment = repelem(cnames, length(t));

df = table(compartment, time);
for i=1:nm
    molecule = mnames{i};
    temp = [];
    for j=1:nc
        temp = [temp; y(:,(j-1)*nm + i)];
    end
    df.(molecule) = temp;
end
writetable(df, "simulation_results.csv");

