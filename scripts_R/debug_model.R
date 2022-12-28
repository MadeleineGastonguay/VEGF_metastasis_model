setwd("~/proj/VEGF_metastasis_model/")
library(tidyverse)
library(here)
library(patchwork)
theme_set(theme_classic())
##### 
# Helper functions
number_per_cell <- function(conc, tissue){
  ESAV = ifelse(tissue == "main", 73, 105)
  conc*6.0230e+23*1e-5/ESAV;
}

convert_to_pmol <- function(conc, tissue){
  K_AV = ifelse(tissue == "primary", 0.58, 0.03)
  conc/K_AV*1e15
}

# Plotting functions
plot_all_receptors <- function(tissue){
  sims %>% filter(compartment == tissue) %>% 
    select(time,  contains("R1"), contains("R2"), contains("N1"), contains("N2")) %>% 
    pivot_longer(!time) %>% 
    mutate(value = number_per_cell(value, tissue)) %>% 
    ggplot(aes(time, value)) + 
    geom_line() + 
    facet_wrap(~name, scales = "free_y") + 
    labs(x = "time (s)", y = "#/cell", title = "All receptor complexes")
}

plot_free_receptors <- function(tissue){
  targets <- data.frame(name = c("R1", "R2", "N1", "N2"), value = c(3750, 300, 20000, 20000))
  
  sims %>% filter(compartment == tissue) %>% 
    select(time, R1, R2, N1, N2) %>% 
    mutate(number_per_cell(across(!time), tissue)) %>%
    pivot_longer(!time) %>% 
    ggplot(aes(time, value)) + 
    geom_line() + 
    facet_wrap(~name, scales = "free_y") + 
    geom_hline(data = targets, aes(yintercept = value), color = "Red") + 
    labs(x = "time (s)", y = "#/cell", title = "Free receptors")
}

plot_non_receptors <- function(tissue){
  K_AV = ifelse(tissue == "primary", 0.58, 0.03)
  
  sims %>% filter(compartment == tissue) %>% 
    select(time,  Mecm, Mebm, Mpbm, Mecm_V165, Mebm_V165, Mpbm_V165, V165) %>% 
    pivot_longer(!time) %>% 
    mutate(value = value/K_AV*1e15, 
           name = factor(name, levels = c("Mecm", "Mebm", "Mpbm", "Mecm_V165", "Mebm_V165", "Mpbm_V165", "V165"))) %>% 
    ggplot(aes(time, value)) + 
    geom_line() + 
    facet_wrap(~name, scales = "free_y") + 
    labs(x = "time (s)", y = "concentration (pM)", title = "Free and matrix-bound ligand")
}

plot_summed_receptors <- function(tissue){
  targets <- data.frame(name = c("R1", "R2", "N1", "N2"), value = c(3750, 300, 20000, 20000))
  
  sims %>% filter(compartment == tissue) %>% 
    select(time, contains("R1"), contains("R2"), contains("N1"), contains("N2")) %>%  
    rowwise %>% 
    mutate(R1 = sum(across(contains("R1"))),
           R2 = sum(across(contains("R2"))),
           N1 = sum(across(contains("N1"))),
           N2 = sum(across(contains("N2")))) %>% 
    mutate(number_per_cell(across(!time), tissue)) %>% 
    select(time, R1, R2, N1, N2) %>% 
    pivot_longer(!time) %>% 
    ggplot(aes(time, value)) + 
    geom_line() + 
    facet_wrap(~name, scales = "free_y") + 
    geom_hline(data = targets, aes(yintercept = value), color = "Red") + 
    labs(x = "time (s)", y = "#/cell", title = "Total number of receptors")
}

quick_bar <- function(df, m,convert = "#/cell", tissue = "primary", position = "stack"){
  converter = switch(convert, "#/cell" = number_per_cell, "pM" = convert_to_pmol)
  y_lab = ifelse(position == "fill", "proportion", convert)
  
  df <- df %>% 
    filter(compartment == tissue, grepl(m, molecule), q_V165 < 1e4) %>% 
    group_by(q_V165) %>% 
    filter(time == max(time)) %>% 
    ungroup
  
  if(m == "V165"){
    df <- df %>% pivot_wider(names_from = molecule, values_from = concentration) %>%  
      group_by(q_V165) %>% 
      summarise(Free = V165, 
                `Bound to M` = sum(across(starts_with("M") & !contains("R"))), 
                MVR = sum(across(starts_with("Mebm_V165_"))),
                RVN = sum(across(starts_with("R2_V165_"))),
                `Receptor Bound` = R1_V165 + R2_V165 + N1_V165 + N2_V165) %>% 
      pivot_longer(!q_V165, names_to = "molecule", values_to = "concentration")
  }
  
  df %>% 
    mutate(value = converter(concentration, tissue)) %>% 
    ggplot(aes("", value, fill = molecule)) + 
    geom_bar(stat = "identity", position = position) + 
    labs(x = "", y = y_lab, fill = "Complex",
         title = str_interp("${m} in ${tissue}")) + 
    facet_wrap(~q_V165, scales = "free_y")
}


#####
# Simulation results with original parameter values
# These results have R2 primarily in the unbound state and very little V165
sims <- read_csv(here("simulation_results_fix_error.csv"))
pdf(here("primary_plots.pdf"))
plot_summed_receptors("primary")
plot_free_receptors("primary")
plot_all_receptors("primary")
plot_non_receptors("primary")

sims %>% filter(compartment == "primary") %>% 
  arrange(time) %>% select(contains("V165")) %>% tail(n = 1) %>% 
  summarise(Free = V165, 
            `Bound to M` = sum(across(starts_with("M") & !contains("R"))), 
            MVR = sum(across(starts_with("Mebm_V165_"))),
            RVN = sum(across(starts_with("R2_V165_"))),
            `Receptor Bound` = R1_V165 + R2_V165 + N1_V165 + N2_V165) %>% 
  pivot_longer(everything()) %>% 
  mutate(value = value/0.58*1e15) %>% 
  ggplot(aes("VEGF165", value, fill = name)) + 
  geom_bar(stat = "identity", position = "stack") + 
  labs(x = "", y = "Concentration (pM)")

sims %>% filter(compartment == "primary") %>% 
  arrange(time) %>% select(contains("R1")) %>% tail(n = 1) %>% 
  pivot_longer(everything()) %>% 
  mutate(value = number_per_cell(value, "primary")) %>% 
  ggplot(aes("R1", value, fill = name)) + 
  geom_bar(stat = "identity", position = "fill") + 
  labs(x = "", y = "#/cell", fill = "Complex") +
  
  sims %>% filter(compartment == "primary") %>% 
  arrange(time) %>% select(contains("R2")) %>% tail(n = 1) %>% 
  pivot_longer(everything()) %>% 
  mutate(value = number_per_cell(value, "primary")) %>% 
  ggplot(aes("R2", value, fill = name)) + 
  geom_bar(stat = "identity", position = "fill") + 
  labs(x = "", y = "#/cell", fill = "Complex") +
  
  sims %>% filter(compartment == "primary") %>% 
  arrange(time) %>% select(contains("N1")) %>% tail(n = 1) %>% 
  pivot_longer(everything()) %>% 
  mutate(value = number_per_cell(value, "primary")) %>% 
  ggplot(aes("N1", value, fill = name)) + 
  geom_bar(stat = "identity", position = "fill") + 
  labs(x = "", y = "#/cell", fill = "Complex") +
  
  sims %>% filter(compartment == "primary") %>% 
  arrange(time) %>% select(contains("N2")) %>% tail(n = 1) %>% 
  pivot_longer(everything()) %>% 
  mutate(value = number_per_cell(value, "primary")) %>% 
  ggplot(aes("N2", value, fill = name)) + 
  geom_bar(stat = "identity", position = "fill") + 
  labs(x = "", y = "#/cell", fill = "Complex")

dev.off()

pdf(here("main_plots.pdf"))
plot_summed_receptors("main")
plot_free_receptors("main")
plot_all_receptors("main")
plot_non_receptors("main")

sims %>% filter(compartment == "main") %>% 
  arrange(time) %>% select(contains("V165")) %>% tail(n = 1) %>% 
  summarise(Free = V165, 
            `Bound to M` = sum(across(starts_with("M") & !contains("R"))), 
            MVR = sum(across(starts_with("Mebm_V165_"))),
            RVN = sum(across(starts_with("R2_V165_"))),
            `Receptor Bound` = R1_V165 + R2_V165 + N1_V165 + N2_V165) %>% 
  pivot_longer(everything()) %>% 
  mutate(value = value*1e12) %>% 
  ggplot(aes("VEGF165", value, fill = name)) + 
  geom_bar(stat = "identity", position = "stack") + 
  labs(x = "", y = "Concentration (pM)")

sims %>% filter(compartment == "main") %>% 
  arrange(time) %>% select(contains("R1")) %>% tail(n = 1) %>% 
  pivot_longer(everything()) %>% 
  mutate(value = number_per_cell(value, "main")) %>% 
  ggplot(aes("R1", value, fill = name)) + 
  geom_bar(stat = "identity", position = "fill") + 
  labs(x = "", y = "#/cell", fill = "Complex") +
  
  sims %>% filter(compartment == "main") %>% 
  arrange(time) %>% select(contains("R2")) %>% tail(n = 1) %>% 
  pivot_longer(everything()) %>% 
  mutate(value = number_per_cell(value, "main")) %>% 
  ggplot(aes("R2", value, fill = name)) + 
  geom_bar(stat = "identity", position = "fill") + 
  labs(x = "", y = "#/cell", fill = "Complex") +
  
  sims %>% filter(compartment == "main") %>% 
  arrange(time) %>% select(contains("N1")) %>% tail(n = 1) %>% 
  pivot_longer(everything()) %>% 
  mutate(value = number_per_cell(value, "main")) %>% 
  ggplot(aes("N1", value, fill = name)) + 
  geom_bar(stat = "identity", position = "fill") + 
  labs(x = "", y = "#/cell", fill = "Complex") +
  
  sims %>% filter(compartment == "main") %>% 
  arrange(time) %>% select(contains("N2")) %>% tail(n = 1) %>% 
  pivot_longer(everything()) %>% 
  mutate(value = number_per_cell(value, "main")) %>% 
  ggplot(aes("N2", value, fill = name)) + 
  geom_bar(stat = "identity", position = "fill") + 
  labs(x = "", y = "#/cell", fill = "Complex")

dev.off()

sims %>% ggplot(aes(time, V165)) + geom_line() + facet_wrap(~compartment)

#####
# Does increasing V165 secretion allow for R2 binding?
# Results show that R2 can bnd when V165 seccretuin rate is 10/cell but that means > 90pM of V165
sec_debug <- read_csv(here("debug_results","simulation_results_debug_secretion.csv"))
str(sec_debug)
sec_debug_long <- sec_debug %>% 
  pivot_longer(!c(q_V165, time), values_to = "concentration",
               names_to = c("compartment", "molecule"), names_sep = "[.]")

pdf(here("primary_at_increasing_qV165.pdf"))
sec_debug %>% 
  filter(q_V165 < 1e4) %>%
  ggplot(aes(time, primary.V165/0.58*1e15, 
             color = as.factor(q_V165  * 1534 / (1e-5 * 6.023e23)))) + 
  geom_line() + 
  labs(x = "time", y = "V165 pM", color = "Secretion (moles/cm^3 tissue)",
       title = "Concentration of V165 in Primary tumor") +
  facet_wrap(~as.factor(q_V165), scales = "free_y") 

quick_bar(sec_debug_long, "V165",  "pM")
quick_bar(sec_debug_long, "R1")
quick_bar(sec_debug_long, "R2")
quick_bar(sec_debug_long, "N1")
quick_bar(sec_debug_long, "N2")

dev.off()

#####
# Does increasing R2 levels induce binding at lower secretion values?
# No!
# archive_increased_R2 <- sec_debug2
sec_debug2 <- read_csv(here("simulation_results_debug_secretion_increasedR2.csv"))
str(sec_debug)
sec_debug_long2 <- sec_debug2 %>% 
  pivot_longer(!c(q_V165, time), values_to = "concentration",
               names_to = c("compartment", "molecule"), names_sep = "[.]")

merge(sec_debug_long, sec_debug_long2, by = c("q_V165", "time", "molecule", "compartment"), 
      suffixes = c("_orig", "_kon")) %>% 
  filter(grepl("R2", molecule), compartment =="primary", q_V165 <= 10) %>%  
  ggplot(aes(number_per_cell(concentration_orig, "primary"), 
             number_per_cell(concentration_kon, "primary"), color = time)) + 
  geom_point() + facet_wrap(~ q_V165 + molecule, scales = "free", ncol = 5) + geom_abline()

pdf(here("primary_at_increasing_qV165_increasedR2.pdf"))
quick_bar(sec_debug_long2, "V165",  "pM") + 
  labs(subtitle = "k_prod.R2 increased by 2 orders of magnitude")
quick_bar(sec_debug_long2, "R1") + 
  labs(subtitle = "k_prod.R2 increased by 2 orders of magnitude")
quick_bar(sec_debug_long2, "R2") + 
  labs(subtitle = "k_prod.R2 increased by 2 orders of magnitude")
quick_bar(sec_debug_long2, "N1") + 
  labs(subtitle = "k_prod.R2 increased by 2 orders of magnitude")
quick_bar(sec_debug_long2, "N2") + 
  labs(subtitle = "k_prod.R2 increased by 2 orders of magnitude")
dev.off()

pdf(here("main_at_increasing_qV165_increasedR2.pdf"))
quick_bar(sec_debug_long2, "V165",  "pM", tissue = "main")
quick_bar(sec_debug_long2, "R1", tissue = "main")
quick_bar(sec_debug_long2, "R2", tissue = "main")
quick_bar(sec_debug_long2, "N1", tissue = "main")
quick_bar(sec_debug_long2, "N2", tissue = "main")
dev.off()



temp_df <- sec_debug_long2 %>% 
  filter(compartment != "blood") %>% 
  pivot_wider(names_from = "molecule", values_from = "concentration") %>% 
  rowwise %>% mutate(R2_total = number_per_cell(sum(across(contains("R2"))), compartment), 
                     R1_total = number_per_cell(sum(across(contains("R1"))), compartment),
                     N1_total = number_per_cell(sum(across(contains("N1"))), compartment), 
                     N2_total = number_per_cell(sum(across(contains("N2"))), compartment)) 

temp_df %>% filter(compartment != "blood") %>% 
  ggplot(aes(time, N2_total, color = as.factor(q_V165))) + 
  geom_line() +
  facet_wrap(~compartment)

temp_df %>% 
  ggplot(aes(time, R1_total, color = as.factor(q_V165))) + 
  geom_line() +
  facet_wrap(~compartment)

## eliminating internalization
sec_debug3 <- read_csv(here("simulation_results_debug_secretion_increasedR2_noint.csv"))
sec_debug_long3 <- sec_debug3 %>% 
  pivot_longer(!c(q_V165, time), values_to = "concentration",
               names_to = c("compartment", "molecule"), names_sep = "[.]")
temp_df3 <- sec_debug_long3 %>% 
  filter(compartment != "blood") %>% 
  pivot_wider(names_from = "molecule", values_from = "concentration") %>% 
  rowwise %>% mutate(R2_total = number_per_cell(sum(across(contains("R2"))), compartment), 
                     R1_total = number_per_cell(sum(across(contains("R1"))), compartment),
                     N1_total = number_per_cell(sum(across(contains("N1"))), compartment), 
                     N2_total = number_per_cell(sum(across(contains("N2"))), compartment)) 

png("only_R2int.png")
temp_df3 %>% filter(compartment != "blood") %>% 
  ggplot(aes(time, R2_total, color = as.factor(q_V165))) + 
  geom_line() +
  facet_wrap(~compartment)
dev.off()

## Calculate internalization values for R2 moleucles
sec_debug %>% 
  filter(q_V165 < 1e3) %>% 
  # filter(compartment == "primary", grepl("R2", molecule), q_V165 < 1e3) %>% 
  ggplot(aes(time, 3.12e-2*primary.R2_V165_N1, color = as.factor(q_V165))) + 
  geom_line()

sec_debug %>% 
  filter(q_V165 < 1e3) %>% 
  mutate(prod_R2V165 = 1e7*1.7241e+03*primary.R1*primary.V165 - 1e-3*primary.R2_V165 + 
           3e13/105*primary.R2_V165*primary.N1 - 1e-3*primary.R2_V165_N1 +
           3e13/105*primary.R2_V165*primary.N2 - 1e-3*primary.R2_V165_N2,
         int_R2V165 = 3.12e-2*primary.R2_V165) %>% 
  # ggplot(aes(prod_R2V165, int_R2V165, color = as.factor(q_V165))) + 
  ggplot(aes(color = as.factor(q_V165))) +
  geom_line(aes(time, prod_R2V165)) + 
  geom_line(aes(time, int_R2V165), lty = "dashed") 


## Simulations with increased R_V internalization
sec_debug4 <- read_csv(here("simulation_results_debug_secretion_increase_RV_int.csv"))
str(sec_debug)
sec_debug_long4 <- sec_debug4 %>% 
  pivot_longer(!c(q_V165, time), values_to = "concentration",
               names_to = c("compartment", "molecule"), names_sep = "[.]")

sec_debug_long4 %>% 
  filter(compartment == "primary", grepl("R2", molecule), q_V165 < 1e4) %>% 
  group_by(q_V165) %>% 
  filter(time == max(time)) %>% 
  ungroup %>% 
  mutate(value = number_per_cell(concentration, "primary")) %>% 
  ggplot(aes("R2", value, fill = molecule)) + 
  geom_bar(stat = "identity", position = "stack") + 
  labs(x = "", y = "#/cell", fill = "Complex") + 
  facet_wrap(~q_V165)

#####
# Does increasing kon R2-V165 induce binding at lower secretion values?
sec_debug5 <- read_csv(here("simulation_results_debug_secretion_increasedKon.csv"))
sec_debug_long5 <- sec_debug5 %>% 
  pivot_longer(!c(q_V165, time), values_to = "concentration",
               names_to = c("compartment", "molecule"), names_sep = "[.]")

pdf(here("primary_at_increasing_qV165_increasedKon.pdf"))
quick_bar(sec_debug_long5, "V165",  "pM")
quick_bar(sec_debug_long5, "R1")
quick_bar(sec_debug_long5, "R2")
quick_bar(sec_debug_long5, "N1")
quick_bar(sec_debug_long5, "N2")
dev.off()

####
# Does removing nerupolin make a difference?
sec_debug6 <- read_csv(here("debug_results","simulation_results_debug_secretion_no_N.csv"))
sec_debug_long6 <- sec_debug6 %>% 
  pivot_longer(!c(q_V165, time), values_to = "concentration",
               names_to = c("compartment", "molecule"), names_sep = "[.]")

pdf(here("fig","primary_at_increasing_qV165_noN.pdf"))
sec_debug6 %>% 
  ggplot(aes(time, primary.V165/0.58*1e15, 
             color = as.factor(q_V165  * 1534 / (1e-5 * 6.023e23)))) + 
  geom_line() + 
  labs(x = "time", y = "V165 pM", color = "Secretion (moles/cm^3 tissue)",
       title = "Concentration of Free V165 in Primary tumor") +
  facet_wrap(~as.factor(q_V165), scales = "free_y") 

quick_bar(sec_debug_long6, "V165",  "pM")
quick_bar(sec_debug_long6, "R1")
quick_bar(sec_debug_long6, "R2")
quick_bar(sec_debug_long6, "N1")
quick_bar(sec_debug_long6, "N2")
dev.off()

comp_df <- merge(
  sec_debug %>% 
    group_by(q_V165) %>% 
    filter(time == max(time)) %>% 
    select(-time) %>% 
    pivot_longer(!q_V165, names_to = c("tissue", "molecule"), names_sep = "[.]",
                 values_to = "orig"),
  
  sec_debug6 %>% 
    group_by(q_V165) %>% 
    filter(time == max(time)) %>% 
    select(-time) %>% 
    pivot_longer(!q_V165, names_to = c("tissue", "molecule"), names_sep = "[.]",
                 values_to = "noN"),
  
  by = c("q_V165", "tissue", "molecule")
)

comp_df %>% 
  filter(tissue == "primary", grepl("V165", molecule)) %>% 
  mutate(orig = 0.58*1e15*orig, noN = 0.58*1e15*noN) %>% 
  ggplot(aes(orig, noN, color = molecule)) + 
  geom_abline(lty = "dashed") + 
  geom_text(aes(label = molecule), show.legend = F) + 
  facet_wrap(~q_V165, scales = "free", labeller = "label_both") + 
  labs(x = "Concentration with Original Model Parameters (pM)",
       y = "Concentration with no Neuropilin (pM)",
       title = "Location of V165")

comp_df %>% 
  filter(tissue == "primary", grepl("R2", molecule), molecule !="R2") %>% 
  mutate(orig = number_per_cell(orig, tissue), noN = number_per_cell(noN, tissue)) %>% 
  ggplot(aes(orig, noN, color = molecule)) + 
  geom_abline(lty = "dashed") + 
  geom_text(aes(label = molecule), show.legend = F) + 
  facet_wrap(~q_V165, scales = "free", labeller = "label_both") + 
  labs(x = "#/cell with Original Model Parameters",
       y = "#/cell with no Neuropilin",
       title = "Location of R2")

comp_df %>% 
  filter(tissue == "primary", grepl("R1", molecule)) %>% 
  mutate(orig = number_per_cell(orig, tissue), noN = number_per_cell(noN, tissue)) %>% 
  ggplot(aes(orig, noN, color = molecule)) + 
  geom_abline(lty = "dashed") + 
  geom_text(aes(label = molecule), show.legend = F) + 
  facet_wrap(~q_V165, scales = "free", labeller = "label_both") + 
  labs(x = "#/cell with Original Model Parameters",
       y = "#/cell with no Neuropilin",
       title = "Location of R1")





