# VEGF_metastasis_model
This repository contains work from my rotation in the Mac Gabhann lab at Johns Hopkins University. We are building a QSP model for VEGF signaling in prostate cancer metastasis when primary and metastatic tumors exhibit different signaling characteristics.

# Repository Contents

**scripts_matlab** contains code for model simulation in addition to code for debuging the model and exploring the effect of varying parameter values. The main scripts for this are multi_tissue_driver_VEGF.m and multi_tissue_VEGF_debug.m, respectively. In addition, this folder contains pseudo code for the stochastic appearance of metastatsis, found in pseudocode_stochastic_simulations.m. 

**scripts_R** contains code and html files from exploring model results in addition to code to generate figures for my lab meeting presentation.

**results** contains concentration-time profiles for each compartment from simulations with updated model paramters.

**debug_results** contains results from simulations with altered parameters. The end of the file name matches a case in scripts_matlab/alter_params.m, indicating which parameter changes were made.

**fig** contains some exploratory figures of model results.

**lab_meeting.pptx** is my presentation from the last lab meeting of the fall 2022 semester. It includes results from the VEGF model in addition to slides exploring the difference between full PBPK models, lumped PBPK models, and compartmental PK models.

**MacGabhann_Rotation_Summary.pdf** is a summary of my work in addition to an outline for future directions of the project.
