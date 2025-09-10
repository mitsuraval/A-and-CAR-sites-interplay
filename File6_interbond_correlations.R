
###########
#corelations

# ─── 0. INSTALL & LOAD Hmisc ────────────────────────────────────────────────────
if (!requireNamespace("Hmisc", quietly = TRUE)) install.packages("Hmisc")
library(Hmisc)
library(tidyverse)

# ─── 1. BUILD YOUR WIDE MATRIX ───────────────────────────────────────────────────
wide_mat <- pca_input %>% 
  select(-structure_id) %>%
  as.matrix()

# ─── 2. GET CORRELATIONS + P-VALUES ─────────────────────────────────────────────
rc    <- rcorr(wide_mat, type = "pearson")
r_mat <- rc$r
p_mat <- rc$P

# ─── 3. MELT INTO A TIDY DF & KEEP UPPER TRIANGLE ───────────────────────────────
corr_df <- as_tibble(r_mat, rownames = "bond1") %>%
  pivot_longer(-bond1, names_to = "bond2", values_to = "r") %>%
  mutate(
    i       = match(bond1, colnames(r_mat)),
    j       = match(bond2, colnames(r_mat)),
    p_value = p_mat[cbind(i, j)]
  ) %>%
  filter(i < j) %>%
  arrange(desc(abs(r)))

# ─── 4. FORMAT p_value WITH 12 SIGNIFICANT DIGITS ────────────────────────────────
corr_df <- corr_df %>%
  mutate(
    p_value = formatC(
      p_value,
      format = "e",   # scientific notation
      digits = 12     # twelve significant digits
    )
  )

# ─── 5. VIEW THE TOP 15 PAIRS ────────────────────────────────────────────────────
print(head(corr_df, 15))

# after you have your corr_df as above...

# 1) Save to disk as CSV
readr::write_csv(corr_df, "bond_correlation_results.csv")

# 2) (If you’re in RStudio) pop it up in the Viewer
View(corr_df)

library(dplyr)

# 1) Define bond groups
A_site    <- c("a","b","f","g","k","l","m","n","o","p","q")
CAR_site  <- c("d","e","i","j","t","u","v","w","x","y","z")

# 2) Filter for cross-site correlations AND p_value < 0.05
cross_corr_df <- corr_df %>%
  # if p_value is character, convert to numeric
  mutate(p_value_num = as.numeric(p_value)) %>%
  filter(
    ((bond1 %in% A_site & bond2 %in% CAR_site) |
       (bond2 %in% A_site & bond1 %in% CAR_site)) &
      p_value_num < 0.05
  ) %>%
  # drop the helper column if you like
  select(-p_value_num)

# 3) Inspect the results
print(cross_corr_df, n = 25)


# ─── 6. PIVOT FOR HEATMAP ─────────────────────────────────────────────────────────
# first, tag each pair so we know which is A_site vs CAR_site
cross_heatmap <- cross_corr_df %>%
  mutate(
    bondA = if_else(bond1 %in% A_site, bond1, bond2),
    bondC = if_else(bond1 %in% CAR_site, bond1, bond2),
    rnum  = as.numeric(r)
  ) %>%
  select(bondA, bondC, rnum)

# ensure the proper ordering on each axis
cross_heatmap <- cross_heatmap %>%
  mutate(
    bondA = factor(bondA, levels = A_site),
    bondC = factor(bondC, levels = CAR_site)
  )

# ─── 7. PLOT HEATMAP ──────────────────────────────────────────────────────────────
library(ggplot2)
ggplot(cross_heatmap, aes(x = bondC, y = bondA, fill = rnum)) +
  geom_tile(color = "grey90") +
  scale_fill_gradient2(
    low      = "steelblue",
    mid      = "white",
    high     = "darkred",
    midpoint = 0,
    limits   = c(-1, 1),
    name     = expression(rho)
  ) +
  labs(
    x     = "CAR-site bond",
    y     = "A-site bond",
    title = "Significant Cross-site Pearson Correlations"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x  = element_text(angle = 45, hjust = 1),
    panel.grid   = element_blank(),
    legend.title = element_text(size = 12),
    legend.text  = element_text(size = 10)
  )

# ─── 8. (OPTIONAL) ANNOTATE WITH r-VALUES ─────────────────────────────────────────
# if you want to print the numeric r inside each tile:
ggplot(cross_heatmap, aes(x = bondC, y = bondA, fill = rnum)) +
  geom_tile(color = "grey90") +
  geom_text(aes(label = sprintf("%.2f", rnum)), size = 3, color = "black") +
  scale_fill_gradient2(
    low      = "steelblue",
    mid      = "white",
    high     = "darkred",
    midpoint = 0,
    limits   = c(-1, 1),
    name     = expression(rho)
  ) +
  labs(
    x     = "CAR-site bond",
    y     = "A-site bond",
    title = "Significant Cross-site Pearson Correlations"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x  = element_text(angle = 45, hjust = 1),
    panel.grid   = element_blank(),
    legend.title = element_text(size = 12),
    legend.text  = element_text(size = 10)
  )

########

# ─── ASSUMING YOU’VE ALREADY RUN STEPS 1–5 AND HAVE corr_df ────────────────
# corr_df has columns bond1, bond2, r (string), p_value (string), plus the helpers i,j

# 1) Define bond groups again
A_site   <- c("a","b","f","g","k","l","m","n","o","p","q")
CAR_site <- c("d","e","i","j","t","u","v","w","x","z")

# 2) Select all cross-site pairs (no p-value filter)
cross_all <- corr_df %>%
  mutate(rnum = as.numeric(r)) %>%
  filter(
    (bond1 %in% A_site & bond2 %in% CAR_site) |
      (bond2 %in% A_site & bond1 %in% CAR_site)
  ) %>%
  # identify which is the A-site vs the CAR-site bond
  mutate(
    bondA = if_else(bond1 %in% A_site, bond1, bond2),
    bondC = if_else(bond1 %in% CAR_site, bond1, bond2)
  ) %>%
  select(bondA, bondC, rnum)

# 3) Ensure proper factor ordering
cross_all <- cross_all %>%
  mutate(
    bondA = factor(bondA, levels = A_site),
    bondC = factor(bondC, levels = CAR_site)
  )

# ─── 4. PLOT THE FULL CROSS-SITE HEATMAP ────────────────────────────────────
library(ggplot2)
ggplot(cross_all, aes(x = bondC, y = bondA, fill = rnum)) +
  geom_tile(color = "grey90") +
  geom_text(
    aes(label = sprintf("%.2f", rnum)),
    size  = 3,
    color = "black"
  ) +
  scale_fill_gradient2(
    low      = "steelblue",
    mid      = "white",
    high     = "darkred",
    midpoint = 0,
    limits   = c(-1,1),
    name     = expression(rho)
  ) +
  labs(
    x     = "CAR-site bond",
    y     = "A-site bond",
    title = "All Cross-site Pearson Correlations"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x  = element_text(angle = 45, hjust = 1),
    panel.grid   = element_blank(),
    legend.title = element_text(size = 12),
    legend.text  = element_text(size = 10)
  )


##########

# ─── 0. INSTALL & LOAD LIBRARIES ─────────────────────────────────────────
if (!requireNamespace("Hmisc", quietly = TRUE)) install.packages("Hmisc")
library(Hmisc)
library(tidyverse)

# ─── 1. BUILD YOUR WIDE MATRIX ─────────────────────────────────────────────
wide_mat <- pca_input %>%
  select(-structure_id) %>%    # drop metadata
  as.matrix()

# ─── 2. GET CORRELATIONS ───────────────────────────────────────────────────
rc    <- rcorr(wide_mat, type = "pearson")
r_mat <- rc$r    # Pearson r
# p_mat <- rc$P  # if you also want p‐values later

# ─── 3. MELT THE FULL MATRIX INTO LONG FORMAT ──────────────────────────────
corr_all <- as_tibble(r_mat, rownames = "bond1") %>%
  pivot_longer(
    -bond1,
    names_to  = "bond2",
    values_to = "r"
  )

# ─── 4. PRESERVE BOND ORDER (OPTIONAL) ─────────────────────────────────────
bond_levels <- rownames(r_mat)
corr_all <- corr_all %>%
  mutate(
    bond1 = factor(bond1, levels = bond_levels),
    bond2 = factor(bond2, levels = rev(bond_levels))  # reverse so (1,1) is top‐left
  )

# ─── 5. PLOT THE HEATMAP WITH OVERLAID r VALUES ────────────────────────────
ggplot(corr_all, aes(x = bond2, y = bond1, fill = r)) +
  geom_tile(color = "grey90") +
  geom_text(aes(label = sprintf("%.2f", r)), size = 3) +
  scale_fill_gradient2(
    low      = "steelblue",
    mid      = "white",
    high     = "darkred",
    midpoint = 0,
    limits   = c(-1, 1),
    name     = "Pearson r"
  ) +
  labs(
    title = "Correlation Matrix of All Bonds",
    x     = "Bond 2",
    y     = "Bond 1"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid  = element_blank()
  )


#######
