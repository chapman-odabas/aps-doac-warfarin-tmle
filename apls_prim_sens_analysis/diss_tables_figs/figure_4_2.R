library(readr)
library(dplyr)
library(ggplot2)

dat <- read_csv("apls_prim_sens_analysis/input_data/tmle_ready/tmle_primary_1y_final_dataset.csv") %>%
  mutate(
    A = as.integer(A),
    treatment_group = factor(A, levels = c(0, 1), labels = c("Warfarin", "DOAC")),
    gender = factor(gender),
    race_clean = factor(race_clean),
    ethnicity_clean = factor(ethnicity_clean)
  )

# If the exported final dataset already excludes these patients, this line is unnecessary.
if ("death_before_index_flag" %in% names(dat)) {
  dat <- dat %>% filter(death_before_index_flag == 0)
}

ps_model <- glm(
  A ~ age_at_index + gender + race_clean + ethnicity_clean +
    htn + dm + dyslipidemia + obesity + tobacco +
    stroke_tia + mi_cad + pad + prior_vte + ckd + liver_dz +
    cancer + copd + heart_failure + atrial_fib + valvular + sle +
    prior_major_bleed + pud + anemia + b2gpi_pos + la_pos +
    arterial_thrombosis + venous_thrombosis + high_risk_phenotype,
  data = dat,
  family = binomial()
)

dat <- dat %>%
  mutate(ps = predict(ps_model, type = "response"))

p_ps <- ggplot(dat, aes(x = ps, fill = treatment_group)) +
  geom_density(alpha = 0.4) +
  labs(
    title = "Propensity score distributions by treatment group",
    x = "Estimated propensity score for DOAC initiation",
    y = "Density",
    fill = "Treatment",
    color = "Treatment"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "right"
  )

print(p_ps)

ggsave(
  filename = "apls_prim_sens_analysis/output_data/final_dissertation_data/tables_figures/figure_4_ps_density.png",
  plot = p_ps,
  width = 9,
  height = 6,
  dpi = 300
)
