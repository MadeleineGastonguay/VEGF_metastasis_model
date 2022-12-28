function [c,m,p,cnames,mnames,pnames,initconc] = declareParams_multi_tissue_VEGF(n_met);
%% Declare tissue compartments (and order in arrays)
c.blood = 1; % Model as written requires the blood compartment to be 1
c.main = 2;
c.primary = 3;
% Add as many metastases as you'd like: 
% (Note that you will have to adjust tissue-specific parameter values
% below or they will automatically be set to the primary tumor value)
for i=1:n_met
    c.(strcat('met', num2str(i))) = i + 3;
end

cnames = fieldnames(c);

%% Declare molecules (and order in arrays)
m.V165 = 1;

m.Mecm = 2;
m.Mebm = 3;
m.Mpbm = 4;

m.R1 = 5;
m.R2 = 6;
m.N1 = 7;
m.N2 = 8;

m.Mecm_V165 = 9;
m.Mebm_V165 = 10;
m.Mpbm_V165 = 11;

m.R1_V165 = 12;
m.R2_V165 = 13;
m.N1_V165 = 14;
m.N2_V165 = 15;

m.R2_V165_N1 = 16;
m.R2_V165_N2 = 17;

m.R1_N1 = 18;
m.R1_N2 = 19;

m.Mebm_V165_R1 = 20;
m.Mebm_V165_R2 = 21;


mnames = fieldnames(m);

%% Establishing Rate Parameters

% BINDING
        
% Ligand:Receptor (source is Lindsey's 2017 paper)
% units in M-1s-1 unless noted otherwise
p.k_on.R1_V165 = 3.0E7;
p.k_on.R2_V165 = 1.0E7;
p.k_on.N1_V165 = 3.125E6; % Bender 2015
p.k_on.N2_V165 = 1E6; % Bender 2015
p.k_on.M_V165 = 1.6E5; %ligand binding to matrix
p.k_on.M_V165_R1 = 3.0E7; %M-bound ligand binding to R1
p.k_on.M_V165_R2 = 1.0E7; %M-bound ligand binding to R2
p.k_on.R1_V165_M = 1.6E5; %M binding to ligand-R1 complex
p.k_on.R2_V165_M = 1.6E5; %M binding to ligand-R2 complex
p.k_on.N1_V165_R2 = 1.0E14; %R2 binding to ligand-N1, unit: (moles/cm^2)-1s-1
p.k_on.R2_V165_N1 = 3.1E13; %N1 binding to ligand-R2, unit: (moles/cm^2)-1s-1
p.k_on.N2_V165_R2 = 1.0E14; %R2 binding to ligand-N1, unit: (moles/cm^2)-1s-1. Assumed same as N1.
p.k_on.R2_V165_N2 = 3.1E13; %N1 binding to ligand-R2, unit: (moles/cm^2)-1s-1. Assumed same as N1.

% units in s-1
p.k_off.R1_V165 = 1.0E-3;
p.k_off.R2_V165 = 1.0E-3;
p.k_off.N1_V165 = 1.0E-3; % Bender 2015
p.k_off.N2_V165 = 1.0E-3; % Bender 2015
p.k_off.M_V165 = 1.0E-2; %ligand binding to matrix
p.k_off.M_V165_R1 = 1.0E-3; %M-bound ligand binding to R1
p.k_off.M_V165_R2 = 1.0E-3; %M-bound ligand binding to R2
p.k_off.R1_V165_M = 1.0E-2; %M binding to ligand-R1 complex
p.k_off.R2_V165_M = 1.0E-2; %M binding to ligand-R2 complex
p.k_off.N1_V165_R2 = 1.0E-3; %R2 binding to ligand-N1
p.k_off.R2_V165_N1 = 1.0E-3; %N1 binding to ligand-R2
p.k_off.N2_V165_R2 = 1.0E-3; %R2 binding to ligand-N1. Assumed same as N1.
p.k_off.R2_V165_N2 = 1.0E-3; %N1 binding to ligand-R2. Assumed same as N1.

% Receptor Coupling
% units: (moles/cm^2 s)^-1
p.k_c.N1_R1 = 1.0E14; 
p.k_c.N2_R1 = 1.0E14; % Bender 2015

% units: (moles/cm^2)-1s-1
p.k_off.N1_R1 = 1.0E-2;
p.k_off.N2_R1 = 1.0E-2; % Bender 2015

% TRAFFICKING and Clearance Parameters 
% Receptor secretion, Endocytosis, Recycling, Degredation
% units: s^-1
p.k_int.R1 = 2.6E-3;
p.k_prod.R1 = 1.7E-16; 

p.k_int.R2 = 2.6E-3;
p.k_prod.R2 = 1.4E-17;

p.k_int.N1 = 2.6E-3;
p.k_prod.N1 = 0.9E-15;

p.k_int.N2 = 2.6E-3; %assume same as N1
p.k_prod.N2 = 0.9E-15;

p.k_int.R2_V = 3.12E-2;
p.k_int.M_V_R2 = 0;
% p.k_int.V_N1_R1 = 3.12E-2; %VEGF165 doesn't bind to N1R1 but 121 does
% p.k_int.V_N2_R1 = 3.12E-2; %VEGF165 doesn't bind to N2R1 but 121 does
p.k_int.M_V_R1 = 0;
p.k_int.R1_V = 3.12E-2;
p.k_int.N1_V = 2.6E-3;
p.k_int.N1_R2_V = 3.12E-2;
p.k_int.N2_V = 2.6E-3; %assume same as N1
p.k_int.N2_R2_V = 3.12E-2; %assume same as N1
p.k_int.N1_R1 = 2.6E-3;
p.k_int.N2_R1 = 2.6E-3; %assume same as N1

% secretion
p.q.V165 = ones(length(cnames), 1)*0.027; % source: Bender 2015, units: #/cell/sec

% Lymphatic drainage: tissue specific 
p.k_l = ones(length(cnames),1)*1.79e-5; %cm^3/s %Source: Feilim claculated it
p.k_l(c.blood) = 0;
p.k_l(c.main) = 0.1418; %cm^3/s %%NOTE: this id from Lindsey's paper, but I think it needs to be re-calculated as we did for the primary tumor

% Permeability parameters: depend on tissue but assumed to be the same for
% each molecule. unit: cm/s
p.k_p = ones(length(cnames),1)*4E-7; %Source: Stefanini 2010
p.k_p(c.blood) = 0;
p.k_p(c.main) = 4.39E-8;

% Clearance from blood: assumed same for all molecules
p.k_cl = 1.08E-3;

%% Geometric Parameters: tissue specific

% Comparment volumes
% units: cm^3
p.U = ones(length(cnames),1)*6.4; %taken from prostate tumor paper
p.U(c.blood) = 5000; %5 Litres = 5000 cm^3
p.U(c.main) = 60453;

% Fraction of interstitial space or plasma
p.K_AV = ones(length(cnames),1)*0.58;%From Bender 2015
p.K_AV(c.blood) = 0.6;
p.K_AV(c.main) = 0.0816; %From Stephanini 2008

% Parameters to convert between measured receptor levels of number per cell
% to moles/cm3 tissue:

% Surface area to volume ratio of endothelial cells
p.ESAV = ones(length(cnames),1)*105; %cm^2/cm^3 tissue % Source: Bender 2015
p.ESAV(c.main) = 73; %cm^2/cm^3 tissue 

% Surface area to volume ratio of parenchyma
p.PSAV = ones(length(cnames),1)*1534; %From Bender 2015
p.PSAV(c.main) = 611;

% Endothelial cell surface area
p.ECSA = ones(length(cnames), 1)*1E-5; %cm^2/EC 

% Avogadro's number
p.N_av = 6.023E23; %molecules/mole

% total abluminal EC surface area: tissue specific
p.S = p.ESAV.*p.U;

% Endothelial cell surface recruitment factor
p.gamma = 1;

% Fraction of EBM accessible to cell surface receptors (tissue specific)
% Assume Mebm wihtin 25nm of cell surface is accessible to surface
% resceptors
% calculate as 25nm/EBM Thickness
p.f = ones(length(cnames),1)*25/50;%EBM thickness taken from Mac Gabhann 2006
p.f(c.main) = 0.286;

pnames = {};
pcat = fieldnames(p);
for i = 1:length(pcat)
    name = pcat{i};
    if isstruct(p.(name))
        for param = fieldnames(p.(name))
            pnames = [pnames; strcat(name, "_", param)];
        end
    else
        pnames = [pnames; name];
    end
end

%% Initial Conditions
for i=1:length(cnames)
    for j=1:length(mnames)
	    initconc.(cnames{i}).(mnames{j}) = 0.;
    end
end

initconc.main.Mecm = 2.15E-11;
initconc.main.Mebm = 2E-12;
initconc.main.Mpbm = 2E-11;
% Need to specify for primary!
initconc.primary.Mecm = 2.15E-11;
initconc.primary.Mebm = 2E-12;
initconc.primary.Mpbm = 2E-11;



return
        
