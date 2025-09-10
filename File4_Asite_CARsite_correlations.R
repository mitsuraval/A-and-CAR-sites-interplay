

###########################################################
###########################################

######CAR site A site hydrogen bonds correlation#######

# Perform both correlation tests
cgu_test <- cor.test(
  combined_data %>% filter(car_site == "CGU") %>% pull(hydrogen_bonds_A),
  combined_data %>% filter(car_site == "CGU") %>% pull(hydrogen_bonds_CAR)
)

gcu_test <- cor.test(
  combined_data %>% filter(car_site == "GCU") %>% pull(hydrogen_bonds_A),
  combined_data %>% filter(car_site == "GCU") %>% pull(hydrogen_bonds_CAR)
)

# Extract raw p-values and estimates
p_raw <- c(cgu_test$p.value, gcu_test$p.value)
estimates <- c(cgu_test$estimate, gcu_test$estimate)
conf_ints <- list(cgu_test$conf.int, gcu_test$conf.int)

# Apply Bonferroni correction
p_bonferroni <- p.adjust(p_raw, method = "bonferroni")

# Create formatted confidence interval strings
ci_strings <- sapply(conf_ints, function(ci) {
  paste0(round(ci[1], 4), " to ", round(ci[2], 4))
})

# Store the result in a data frame
correlation_results_df <- tibble(
  car_site = c("CGU", "GCU"),
  correlation = round(estimates, 4),
  p_value_raw = signif(p_raw, 3),
  p_value_bonferroni = signif(p_bonferroni, 3),
  `95% CI` = ci_strings,
  significant = p_bonferroni < 0.05
)

print(correlation_results_df)

#visualise the correlation

library(ggplot2)

ggplot(combined_data, aes(x = hydrogen_bonds_A, y = hydrogen_bonds_CAR, color = car_site)) +
  geom_point(alpha = 0.5, size = 3) +  # Larger points
  geom_smooth(method = "lm", se = FALSE, linewidth = 1.2) +  # Thicker lines
  labs(
    x = "A-site Hydrogen Bonds",
    y = "CAR-site Hydrogen Bonds",
    title = "Correlation between A-site and CAR-site Hydrogen Bonds",
    color = "CAR Site"
  ) +
  theme_minimal(base_size = 16) +  # Increase base font size
  theme(
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 18, face = "bold"),
    axis.text = element_text(size = 14),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14)
  )

#######alternate stratified A-site correlations visualization

library(dplyr)
library(tibble)
library(ggplot2)
library(purrr)

# ─── 1. COMPUTE CORRELATIONS FOR BOTH GEOMETRIES ───────────────────────────────
correlation_all <- combined_data %>%
  group_by(w_position, a_site, car_site) %>%
  summarise(
    cor_test = list(cor.test(hydrogen_bonds_A, hydrogen_bonds_CAR)),
    .groups = "drop"
  ) %>%
  mutate(
    correlation    = map_dbl(cor_test, ~ .x$estimate),
    p_value_raw    = map_dbl(cor_test, ~ .x$p.value),
    conf_low       = map_dbl(cor_test, ~ .x$conf.int[1]),
    conf_high      = map_dbl(cor_test, ~ .x$conf.int[2])
  ) %>%
  select(-cor_test) %>%
  mutate(
    p_value_bonferroni = p.adjust(p_value_raw, method = "bonferroni"),
    significant        = p_value_bonferroni < 0.05
  )




################################################
# Stratified correlation analysis and visualization
# Perform correlation within each (A-site codon, CAR-site codon) group



# First group the data
correlation_results_by_a_site <- combined_data %>%
  group_by(a_site, car_site) %>%
  summarise(
    cor_test = list(cor.test(hydrogen_bonds_A, hydrogen_bonds_CAR)),
    .groups = "drop"
  ) %>%
  mutate(
    correlation = map_dbl(cor_test, ~ .x$estimate),
    p_value_raw = map_dbl(cor_test, ~ .x$p.value),
    conf_low = map_dbl(cor_test, ~ .x$conf.int[1]),
    conf_high = map_dbl(cor_test, ~ .x$conf.int[2])
  ) %>%
  select(-cor_test)

# Bonferroni correction across all tests
correlation_results_by_a_site <- correlation_results_by_a_site %>%
  mutate(
    p_value_bonferroni = p.adjust(p_value_raw, method = "bonferroni"),
    significant = p_value_bonferroni < 0.05
  )

# View the results
print(correlation_results_by_a_site)

# Optional: Save to CSV
# write.csv(correlation_results_by_a_site, "correlation_by_a_site.csv", row.names = FALSE)

# -------------------------
# Visualization 
# -------------------------
# Plot correlations by A-site codon and CAR-site codon context
ggplot(correlation_results_by_a_site, aes(x = a_site, y = correlation, color = car_site)) +
  geom_point(size = 4) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(
    title = "Correlation between A-site and CAR-site Hydrogen Bonds by A-site Codon",
    x = "A-site Codon",
    y = "Pearson Correlation (r)",
    color = "+1 Codon Context"
  ) +
  theme_minimal(base_size = 16) +
  theme(
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 18, face = "bold"),
    axis.text = element_text(size = 14),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14)
  )

