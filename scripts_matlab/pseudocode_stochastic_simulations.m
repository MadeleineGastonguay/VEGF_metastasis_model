% This script contains pseudo code for introducing stochastic appearance of
% metastases

%% Start simulation without metastatses
n_met = 0;% Define the number of metastases
[c,m,p,cnames,mnames,pnames,initconc] = declareParams_multi_tissue_VEGF(n_met);
nm = length(mnames);
nc = length(cnames);

% CHANGE: Set t_end to the duration of simulation before using a
% stochastic process to determine if a metastatsis has formed
t_end = 2;

y0=zeros(nm*nc,1);
for i=1:nc
    for j=1:nm
        index = (i-1)*nm + j;
        y0(index)=initconc.(cnames{i}).(mnames{j});
    end
end

[T0,Y0] = multi_tissue_main_VEGF(c,p,m,y0,t_end,runVar);

%% Determine if a metastasis has formed

%
% ADD CODE for stochastic process.
%

%% If it is determined that there is a new metastasis:
n_met = 1;% Define the number of metastases
[c,m,p,cnames,mnames,pnames,initconc] = declareParams_multi_tissue_VEGF(n_met);
nm = length(mnames);
nc = length(cnames);

% CHANGE: Set t_end to the duration of simulation before using a
% stochastic process to determine if a metastatsis has formed
t_end = 2;

y0=zeros(nm*nc,1);
% Set end of last simulations as starting point for these simulations
y0(1:size(Y0,2)) = Y0(end,:);
% Initialize new metastasis compartment
for j=1:nm
    index = (nc-1)*nm + j;
    y0(index)=initconc.(cnames{nc}).(mnames{j});
end

[T1,Y1] = multi_tissue_main_VEGF(c,p,m,y0,t_end,runVar);

%% If it is determined that there is not a new metastasis:
% Don't change any of the parameters declared in the first chunk:

% CHANGE: Set t_end to the duration of simulation before using a
% stochastic process to determine if a metastatsis has formed
t_end = 2;

% Use end of last simulation as initial condition for this simulation
y0 = Y0(end,:);
[T,Y] = multi_tissue_main_VEGF(c,p,m,y0,t_end,runVar);
T0 = [T0; T];
Y0 = [Y0; Y];

%% Determine if a metastasis has formed
% Return to line 24 and repeat. This should be wrapped into a function that
% takes the total simulation duration (simulation_end), times at which metastases may form,
% and the maximum number of metastases. In addition,
% declareParams_multi_tissue_VEGF.m should be updated if metastasis
% parameters vary from primary tumors.

%% Compile results
% Y0, Y1, Y2, ... will all have a different number of columns because they are from
% simulations with different number of compartments. To fix this, we set
% the concentration for metastatses to zero when there is no metastasis:

max_met = 8; % Define the maximum number of metastases
Y0 = [Y0, zeros(size(Y0, 1), nm*(max_met - size(Y0, 2)/nm))];
Y1 = [Y1, zeros(size(Y1, 1), nm*(max_met - size(Y1, 2)/nm))];
%...

% All Ys should now have (3 + max_met)*nm columns and can be combined:
Y = [Y0; Y1; Y2];

% T is the time from start of simulation to end, incremented by the
% interval at which we asked the solver to report concentrations (in this
% case, 60 seconds):
T = 0:60:simulation_end;

%% Notes
% This script outlines one approach. It may be easier to run the simulation
% with the max numer of metastatses at every time step, but "zero-out" the
% metastasis compartments by setting all their parameters equal to zero
% unless the stochastic process indicates the appearance of a new
% metastatsis.

