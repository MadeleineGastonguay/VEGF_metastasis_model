#####
setwd("~/proj/VEGF_metastasis_model/")
library(tidyverse)
library(here)
library(patchwork)
library(ggnewscale)
library(cowplot)
theme_set(theme_classic(base_size = 18))

#####
# Read in the data

# Original parameters
fluxes <- read_csv(here("debug_results", "calculate_fluxes.csv")) %>%
  mutate(q_V165 = factor(q_V165), kprod_R2 = factor(kprod_R2))

orig_simulations <- read_csv(here("debug_results", "simulation_results_debug_R2.csv"))
orig_simulations_long <- orig_simulations %>%
  pivot_longer(!c(q_V165, kprod_R2, time), values_to = "concentration",
               names_to = c("compartment", "molecule"), names_sep = "[.]")

# Updated parameters
new_param_fluxes <- read_csv(here("debug_results", "calculate_fluxes_new_lymph.csv")) %>%
  mutate(q_V165 = factor(q_V165), kprod_R2 = factor(kprod_R2))

new_simulations <- read_csv(here("debug_results", "simulation_results_debug_R2_new_lymph.csv"))
new_simulations_long <- new_simulations %>%
  pivot_longer(!c(q_V165, kprod_R2, time), values_to = "concentration",
               names_to = c("compartment", "molecule"), names_sep = "[.]")

# Simulations with no ligand
simulations_no_ligand <- read_csv(here("debug_results", "simulation_results_debug_R2_no_ligand_new_lymph.csv"))

simulations_no_ligand_long <- simulations_no_ligand %>%
  pivot_longer(!c(q_V165, kprod_R2, time), values_to = "concentration",
               names_to = c("compartment", "molecule"), names_sep = "[.]")
##### 
#Define functions to convert to pM or to #/cell
number_per_cell <- function(conc, tissue){
  ESAV = ifelse(tissue == "main", 73, 105)
  conc*6.0230e+23*1e-5/ESAV;
}

convert_to_pmol <- function(conc, tissue){
  K_AV = ifelse(tissue == "primary", 0.58, 0.03)
  conc/K_AV*1e15
}

convert_secretion_mol_per_cm3_tissue <- function(q_V165, tissue = "primary"){
  q_V165  * 1534 / (1e-5 * 6.023e23)
}

#####
# Plotting functions
quick_bar <- function(df_long, m, tissue = "primary", pos = "stack", convert = "pM", x = NULL, scales = "free_y", outline = NULL){
  if(is.null(x)) x <- ifelse(m == "V165", "kprod_R2", "q_V165")
  facet = c("kprod_R2", "q_V165")[c("kprod_R2", "q_V165") != x]
  
  x_lab = switch(x, "kprod_R2" = "R2 Production (moles/cm^3 tissue)", "q_V165" = "V165 Secretion (#/cell)")
  y_lab = ifelse(pos == "fill", "Fraction", convert)
  # y_lab = switch(pos, "stack" = str_interp("${m} (${convert})"), "fill" = str_interp("Fraction of ${m}"))
  converter = switch(convert, "pM" = convert_to_pmol, "#/cell" = number_per_cell)
  title = switch(tissue, "primary" = str_interp("${m} in Primary Tumor"), 
                 "main" = str_interp("${m} in Main Body"),
                 "blood" = str_interp("${m} in Blood"))
  
  df_long <- df_long %>%
    group_by(q_V165, kprod_R2) %>%
    filter(time == max(time), grepl(m, molecule), compartment == tissue) 
  
  if(m == "V165" & is.null(outline)){
    df_long <- df_long %>% 
      ungroup %>% 
      pivot_wider(names_from = molecule, values_from = concentration) %>%  
      group_by(q_V165, kprod_R2) %>% 
      summarise(Free = V165, 
                `Bound to M` = sum(across(starts_with("M") & !contains("R"))), 
                MVR = sum(across(starts_with("Mebm_V165_"))),
                # RVN = sum(across(starts_with("R2_V165_"))),
                `Receptor Bound` = R1_V165 + R2_V165 + N1_V165 + N2_V165 + R2_V165_N1 + R2_V165_N2) %>% 
      pivot_longer(!c(q_V165, kprod_R2), names_to = "molecule", values_to = "concentration")
  }
  
  p <- df_long %>%
    ggplot(aes_string(str_interp("as.factor(${x})"), "converter(concentration, tissue)", fill = "molecule")) +
    facet_wrap(as.formula(paste0("~", facet)), scales = scales, nrow = 1, labeller = "label_both") +
    labs(x = x_lab, y = y_lab, title = title)
  
  if(!is.null(outline)){
    p + geom_bar(stat = "identity", pos = pos, aes(color = grepl(outline, molecule)))  +
      scale_color_manual(values = c("TRUE" = "black", "FALSE" = NULL)) +
      labs(color = "R2 Complex")
  }else{
    p + geom_bar(stat = "identity", pos = pos)
  }
}

plot_norm_flux <- function(flux){
  flux %>%
    pivot_longer(4:6) %>%
    mutate(norm_value = abs(value/Secretion)) %>%
    ggplot(aes(kprod_R2, norm_value, color = name)) +
    geom_point(size = 3) +
    facet_wrap(~q_V165, labeller= "label_both")  +
    labs(x = "Production of R2 (moles/cm^3 tissue)", y = "Normalized Flux", color = "Flux") +
    theme(legend.position = c(1,0), legend.justification = c(1,0))
}

##### 
# Bar plots of steady state concentrations for original model parameters 
quick_bar(orig_simulations_long %>% filter(q_V165 == 0.027, kprod_R2 == 1.4e-17), "V165", convert = "pM") +
  theme(legend.position = "right", legend.justification = c(0, 0.5), legend.background = element_blank(),
        strip.background = element_blank(), strip.text = element_blank(),
        axis.text.x = element_blank(), axis.title.x = element_blank(), axis.ticks.x = element_blank()) +
  labs(filll = "")

quick_bar(orig_simulations_long %>% filter(q_V165 == 0.027, kprod_R2 == 1.4e-17), "R1", convert = "#/cell") + 
  quick_bar(orig_simulations_long %>% filter(q_V165 == 0.027, kprod_R2 == 1.4e-17), "R2", convert = "#/cell") +
  quick_bar(orig_simulations_long %>% filter(q_V165 == 0.027, kprod_R2 == 1.4e-17), "N1", convert = "#/cell") +
  quick_bar(orig_simulations_long %>% filter(q_V165 == 0.027, kprod_R2 == 1.4e-17), "N2", convert = "#/cell") +
  plot_layout(ncol = 2) &
  theme(legend.position = "right", legend.justification = c(0.5, 0.5), legend.background = element_blank(),
        strip.background = element_blank(), strip.text = element_blank(),
        axis.text.x = element_blank(), axis.title.x = element_blank(), axis.ticks.x = element_blank())

#####
# Bar plots of steady state concentrations for original model parameters with increasing R2 production and VEGF secretion
quick_bar(orig_simulations_long, "R2", convert = "#/cell")
quick_bar(orig_simulations_long %>% filter(molecule != "R2"), "R2", convert = "#/cell")
quick_bar(orig_simulations_long, "V165", convert = "pM")+
  theme(axis.text.x = element_text(angle = 20, hjust = 1))

#####
# Fluxes with original model parameters
fluxes %>% 
  filter(q_V165 == "0.027", kprod_R2 == "1.4e-17") %>% 
  mutate(across(where(is.numeric))/Secretion) %>% 
  pivot_longer(!c(q_V165, kprod_R2)) %>% 
  mutate(name = factor(name, labels = c("Internalization", "Lymphatic Drainage", "Permeability", "Secretion"))) %>% 
  ggplot(aes(name, abs(value))) + 
  geom_bar(stat="identity") +
  theme(axis.text.x = element_text(angle = 20, hjust = 1)) + 
  labs(x = "", y = "Normalized Flux")

# Fluxes with updated lymphatic drainage
new_param_fluxes %>% 
  filter(q_V165 == "0.027", kprod_R2 == "1.4e-17") %>% 
  mutate(across(where(is.numeric))/Secretion) %>% 
  pivot_longer(!c(q_V165, kprod_R2)) %>% 
  mutate(name = factor(name, labels = c("Internalization", "Lymphatic Drainage", "Permeability", "Secretion"))) %>% 
  ggplot(aes(name, abs(value))) + 
  geom_bar(stat="identity") +
  theme(axis.text.x = element_text(angle = 20, hjust = 1)) + 
  labs(x = "", y = "Normalized Flux")

#####
# Bar plots of steady state concentrations for updated model parameters with increasing R2 production and VEGF secretion
quick_bar(new_simulations_long, "R2", convert = "#/cell")
quick_bar(new_simulations_long %>% filter(molecule != "R2"), "R2", convert = "#/cell")
quick_bar(new_simulations_long, "V165", convert = "pM")+
  theme(axis.text.x = element_text(angle = 20, hjust = 1))

#####
# Compare concentrations from simulations with original and updated parameters
categories <- c(
  "Mebm" = "Free Matrix",
  "Mecm" = "Free Matrix",
  "Mpbm" = "Free Matrix",
  "Mebm_V165" = "Matrix-bound VEGF",
  "Mecm_V165" = "Matrix-bound VEGF",
  "Mpbm_V165" = "Matrix-bound VEGF",
  "Mebm_V165_R1" = "MVR",
  "Mebm_V165_R2" = "MVR",
  "N1" = "Free N1",
  "N2" = "Free N2",
  "R1" = "Free R1",
  "R2" = "Free R2",
  "R1_N1" = "Coupled R1-N",
  "R1_N2" = "Coupled R1-N",
  "R1_V165" = "R1-bound VEGF",
  "R2_V165" = "R2-bound VEGF",
  "N1_V165" = "N1-bound VEGF",
  "N2_V165" = "N2-bound VEGF",
  "R2_V165_N1" = "R2-bound VEGF",
  "R2_V165_N2" = "R2-bound VEGF",
  "V165" = "Free VEGF")


compare_df <- merge(
  orig_simulations_long %>% 
    filter(time == max(time), compartment == "primary") %>% 
    mutate(concentration = convert_to_pmol(concentration, "primary")),
  
  new_simulations_long %>% 
    filter(time == max(time), compartment == "primary") %>% 
    mutate(concentration = convert_to_pmol(concentration, "primary")), 
  
  by = names(orig_simulations_long)[1:5],
  suffixes = c("_orig", "_new")
) %>% 
  pivot_longer(starts_with("concentration"), names_to = "simulation", values_to = "concentration")


totals <- compare_df %>% group_by(q_V165, kprod_R2, simulation) %>% 
  summarise(total_V165 = sum(concentration*grepl("V165", molecule)),
            total_R2 = sum(concentration*grepl("R2", molecule)),
            total_R1 = sum(concentration*grepl("R1", molecule)),
            total_N1 = sum(concentration*grepl("N1", molecule)),
            total_N2 = sum(concentration*grepl("N2", molecule))) %>% 
  pivot_longer(starts_with("total"), names_to = "category") %>% 
  pivot_wider(names_from = "simulation")
  

color = "log(as.numeric(as.character(kprod_R2))/convert_secretion_mol_per_cm3_tissue(as.numeric(as.character(q_V165)))))"
p1 <- compare_df %>% 
  mutate(category = categories[molecule]) %>%
  group_by(q_V165, kprod_R2, simulation, category) %>% 
  summarise(concentration = sum(concentration)) %>% 
  pivot_wider(names_from = simulation, values_from = concentration) %>% 
  filter(grepl("bound", category) & !grepl("Matrix", category)) %>%
  ggplot(aes_string("concentration_orig", "concentration_new",
             color = color)) + 
  geom_point(size = 2) + 
  geom_abline() + 
  facet_wrap(~category, scales = "free", nrow = 1) +
  scale_x_log10() + scale_y_log10() +
  scale_color_gradient2() + 
  labs(color = "", x = "", y= "") + 
  theme(axis.text.x = element_text(angle = 25, hjust = 1))

p2 <- totals %>% 
  filter(category != "total_V165") %>%
  mutate(category = gsub("total_", "Total ", category)) %>% 
  ggplot(aes_string("concentration_orig", "concentration_new",
                    color = color)) + 
  geom_point(size = 2) + 
  geom_abline() + 
  facet_wrap(~category, scales = "free", nrow = 1) +
  scale_x_log10() + scale_y_log10() +
  scale_color_gradient2() + 
  labs(color = "", x = "", y= "") + 
  theme(axis.text.x = element_text(angle = 25, hjust = 1))


p3 <- compare_df %>% 
  mutate(category = categories[molecule]) %>%
  group_by(q_V165, kprod_R2, simulation, category) %>% 
  summarise(concentration = sum(concentration)) %>% 
  pivot_wider(names_from = simulation, values_from = concentration) %>% 
  rbind(totals) %>% 
  filter(category %in% c("Matrix-bound VEGF", "MVR", "total_V165")) %>%
  mutate(category = gsub("total_V165", "Total VEGF", category)) %>% 
  ggplot(aes_string("concentration_orig", "concentration_new",
                    color = color)) + 
  geom_point(size = 2) + 
  geom_abline() + 
  facet_wrap(~category, scales = "free", nrow = 1) +
  scale_x_log10() + scale_y_log10() +
  scale_color_gradient2() + 
  labs(color = expression(log(frac("R2 Production","VEGF Secretion"))), x = "", y= "") + 
  theme(axis.text.x = element_text(angle = 25, hjust = 1))

plot_grid(plot_grid(p3 + theme(legend.position = "none"), get_legend(p3), rel_widths = c(3,1)), 
          p1+ theme(legend.position = "none"), p2+ theme(legend.position = "none"), 
          nrow = 3)

#####
# Contour Plots for sensitivity analysis with updated model parameters
contour_df <- new_simulations_long %>% 
  filter(compartment == "primary") %>% 
  group_by(q_V165, kprod_R2) %>% 
  filter(time == max(time)) %>% 
  summarise(total_R2 = number_per_cell(sum(concentration*grepl("R2", molecule)), "primary"),
            free_V165 = convert_to_pmol(sum(concentration*(molecule == "V165")), "primary"),
            bound_R2 = number_per_cell(sum(concentration*(grepl("R2", molecule) & grepl("V165", molecule))), "primary"),
            perc_bound_R2 = bound_R2/total_R2*100)

contour_df %>% 
  # convert secretion to moles/cm^3 tissue
  ggplot(aes(kprod_R2, q_V165  * 1534 / (1e-5 * 6.023e23))) + 
  # Contours for total R2
  geom_contour(aes(z = log(total_R2), color = after_stat(level),
                   # Mark 3000 R2:
                   lty = after_stat(level) == after_stat(level)[which.min(abs(after_stat(level) - log(3000)))]),
               bins = 10 ) +
  scale_color_gradient2("Total R2 (#/cell)", 
                        low = "#762A83", mid = "white", high = "#1B7837",#) +
                        breaks = log(c(3, 30, 3e2, 3e3, 3e4, 3e5)),
                        labels = function(breaks) round(exp(breaks))) +
  # Add contours for free VEGF
  new_scale("color") +
  geom_contour(aes(z = log(free_V165), color = after_stat(level)),
               bins = 10) +
  scale_color_gradient("Free V165 (pM)",
                       low = "#e8cced", high = "#762A83",
                       breaks = log(c(0,1,10,100,1000,10000,100000)),
                       labels = function(breaks) round(exp(breaks))) +
  # Formatting
  labs(x = "Production of R2 in moles/cm^3 tissue/s", y = "V165 Secretion in moles/cm^3 tissue/s",
       title = "Contour plot of Total R2 and Free V165", 
       subtitle = "Log Scale",
       lty = "Total R2 near 3000") +
  scale_x_log10() + scale_y_log10()

# Amount of bound R2
contour_df %>% 
  ggplot(aes(kprod_R2, q_V165  * 1534 / (1e-5 * 6.023e23))) + 
  geom_contour(aes(z = log(bound_R2), color = after_stat(level)),
               bins = 10 ) +
  scale_color_gradient2("Bound R2 (#/cell)", 
                        low = "#762A83", mid = "white", high = "#1B7837",
                        breaks = log(c(3, 30, 3e2, 3e3, 3e4, 3e5)),
                        labels = function(breaks) round(exp(breaks))) +
  new_scale("color") +
  geom_contour(aes(z = log(free_V165), color = after_stat(level)),
               bins = 10) +
  scale_color_gradient("Free V165 (pM)",
                       low = "#e8cced", high = "#762A83",
                       breaks = log(c(0,1,10,100,1000,10000,100000)),
                       labels = function(breaks) round(exp(breaks))) +
  labs(x = "Production of R2 in moles/cm^3 tissue/s", y = "V165 Secretion in moles/cm^3 tissue/s",
       title = "Contour plot of Bound R2 and Free V165", subtitle = "Log Scale", lty = "Total R2 near 3000") +
  scale_x_log10() + scale_y_log10()

# Percent of R2 bound 
contour_df %>% 
  ggplot(aes(kprod_R2, q_V165  * 1534 / (1e-5 * 6.023e23))) + 
  geom_contour(aes(z = perc_bound_R2, color = after_stat(level)),
               bins = 20) +
  scale_color_gradient2("R2 bound to\nLigand (%)",
                        low = "#659473", high = "#1B7837") +
  labs(x = "Production of R2 in moles/cm^3 tissue/s", y = "V165 Secretion in moles/cm^3 tissue/s",
       title = "Contour plot of Percent Bound R2",
       subtitle = "Log Scale") +
  scale_x_log10() + scale_y_log10()


# Overlay R2 with ligand and without
contour_df2 <- new_simulations_long %>% 
  filter(compartment == "primary", time == max(time)) %>% 
  merge(
    simulations_no_ligand_long %>% 
      filter(compartment == "primary"),
    by = c("kprod_R2", "time", "compartment", "molecule")
  ) %>% 
  group_by(q_V165.x, kprod_R2) %>% 
  summarise(total_R2.x = number_per_cell(sum(concentration.x*grepl("R2", molecule)), "primary"),
            free_V165.x = convert_to_pmol(sum(concentration.x*(molecule == "V165")), "primary"),
            perc_bound_R2.x = number_per_cell(sum(concentration.x*(grepl("R2", molecule) & grepl("V165", molecule))), "primary")/total_R2.x*100,
            total_R2.y = number_per_cell(sum(concentration.y*grepl("R2", molecule)), "primary")) %>% 
  mutate(ratio_total_R2 = total_R2.x/total_R2.y)

contour_df2 %>% 
  ggplot(aes(kprod_R2, q_V165.x  * 1534 / (1e-5 * 6.023e23))) + 
  geom_contour(aes(z = log(total_R2.x), color = after_stat(level),
                   # Mark 3000 R2:
                   lty = after_stat(level) == after_stat(level)[which.min(abs(after_stat(level) - log(3000)))]),
               bins = 10 ) +
  scale_color_gradient2("Total R2 with Ligand (#/cell)", 
                        low = "#762A83", mid = "white", high = "#1B7837",#) +
                        breaks = log(c(3, 30, 3e2, 3e3, 3e4, 3e5)),
                        labels = function(breaks) round(exp(breaks))) +
  
  new_scale("color") +
  geom_contour(aes(z = log(total_R2.y), color = after_stat(level),#),
               lty = after_stat(level) == after_stat(level)[which.min(abs(after_stat(level) - log(3000)))]),
               bins = 10) +
  scale_color_gradient("Total R2 without Ligand (#/cell)",
                       low = "#edccd4", high = "#832a39",
                       breaks = log(c(0,1,10,100,1000,10000,100000)),
                       labels = function(breaks) round(exp(breaks))) +
  labs(x = "Production of R2 in moles/cm^3 tissue/s", y = "V165 Secretion in moles/cm^3 tissue/s",
       title = "Contour plot of Total R2 in Simulations With and Without Ligand", 
       subtitle = "Log Scale",
       lty = "Total R2 near 3000") +
  scale_x_log10() + scale_y_log10()

# Contour plot for the ratio of R2 with and without ligand
contour_df2 %>% 
  ggplot(aes(kprod_R2, q_V165.x  * 1534 / (1e-5 * 6.023e23))) + 
  geom_contour(aes(z = ratio_total_R2, color = after_stat(level)),
               bins = 10 ) +
  scale_color_gradient2("Ratio of total R2\nwith ligand vs without", 
                        low = "#762A83", mid = "white", high = "#1B7837") + 
  labs(x = "Production of R2 in moles/cm^3 tissue/s", y = "V165 Secretion in moles/cm^3 tissue/s",
       title = "Contour plot of the Ratio of Total R2 from Model with Ligand to Model Without Ligand", 
       subtitle = "Log Scale") +
  scale_x_log10() + scale_y_log10()
