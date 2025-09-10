
#plot data distribution to test normality
# Fit a preliminary ANOVA model
model <- aov(stack_value ~ car_site * w_position * a_site, data = stack_pair_data)
anova_result <- aov(stack_value ~ car_site * w_position * a_site, data = stack_pair_data)
summary(anova_result)


#########################################################
#normality diagnostics plots for Hbonds

# ─── 0. PACKAGES ───────────────────────────────────────────────────────────────
if (!requireNamespace("tidyverse", quietly=TRUE))  install.packages("tidyverse")
if (!requireNamespace("patchwork", quietly=TRUE))  install.packages("patchwork")
if (!requireNamespace("ggtext", quietly=TRUE))     install.packages("ggtext")
library(tidyverse)
library(patchwork)
library(ggtext)

# ─── 1. FIT YOUR ANOVA MODELS ─────────────────────────────────────────────────
mod_A <- aov(
  hydrogen_bonds ~ car_site * w_position * a_site,
  data = A_anova_long
)
mod_C <- aov(
  hydrogen_bonds ~ car_site * w_position * a_site,
  data = CAR_anova_long
)

# ─── 2. EXTRACT RESIDUALS & RUN SHAPIRO TESTS ─────────────────────────────────
resid_A <- residuals(mod_A)
resid_C <- residuals(mod_C)
sw_A    <- shapiro.test(resid_A)
sw_C    <- shapiro.test(resid_C)

# ─── 3. BUILD A LONG DF FOR PLOTTING ────────────────────────────────────────────
df_res <- tibble(
  residual = c(resid_A, resid_C),
  model    = rep(c("A-site\nH-bonds","CAR-site\nH-bonds"),
                 times = c(length(resid_A), length(resid_C)))
)

sw_df <- tibble(
  model = c("A-site\nH-bonds","CAR-site\nH-bonds"),
  label = c(
    paste0("SW p = ", signif(sw_A$p.value, 2)),
    paste0("SW p = ", signif(sw_C$p.value, 2))
  )
)

# ─── 4. DEFINE FACET TITLE LABELS WITH COLORED TEXT ────────────────────────────
# Match your stack‐figure colors: red for tRNA‐loop & mRNA (N1_N2 etc), blue for CAR
facet_labels <- c(
  "A-site\nH-bonds"  = "<span style='color:#1F78B4;font-size:14pt;'>A-site<br>H-bonds</span>",
  "CAR-site\nH-bonds"= "<span style='color:#1F78B4;font-size:14pt;'>CAR-site<br>H-bonds</span>"
)

# ─── 5. HISTOGRAM + DENSITY PANEL ──────────────────────────────────────────────
p_hist <- ggplot(df_res, aes(residual)) +
  geom_histogram(aes(y = ..density..),
                 bins      = 30,
                 fill      = "#A6CEE3",
                 color     = "#1F78B4",
                 size      = 0.3) +
  geom_density(color = "#1F78B4", size = 0.8) +
  facet_wrap(~model, nrow = 1, labeller = labeller(model = facet_labels)) +
  geom_text(
    data      = sw_df,
    aes(x = Inf, y = Inf, label = label),
    hjust     = 1.1, vjust = 1.2,
    size      = 4, color = "#333333",
    fontface  = "italic"
  ) +
  labs(x = "Residual", y = "Density") +
  theme_minimal(base_size = 14) +
  theme(
    strip.text         = element_markdown(margin = margin(b = 8)),
    panel.grid.major.y = element_line(color = "#eeeeee"),
    panel.grid.minor   = element_blank()
  )

# ─── 6. Q–Q PANEL ───────────────────────────────────────────────────────────────
p_qq <- ggplot(df_res, aes(sample = residual)) +
  stat_qq(color = "#000000", size = 1) +
  stat_qq_line(linetype = "dashed", color = "#E31A1C", size = 0.6) +
  facet_wrap(~model, nrow = 1, labeller = labeller(model = facet_labels)) +
  labs(x = "Theoretical Quantiles", y = "Sample Quantiles") +
  theme_minimal(base_size = 14) +
  theme(
    strip.text         = element_blank(),  # titles already on p_hist
    panel.grid.major.y = element_line(color = "#eeeeee"),
    panel.grid.minor   = element_blank()
  )

# ─── 7. COMBINE & SAVE ──────────────────────────────────────────────────────────
final <- (
  p_hist / p_qq
) + plot_layout(heights = c(1, 1)) +
  plot_annotation(
    title = "Residual Normality Diagnostics for H-bond ANOVAs",
    theme = theme(
      plot.title = element_text(face="bold", size=18, hjust=0.5)
    )
  )

print(final)


##################
#normality diagnostics plots for stacks

library(dplyr)
library(ggplot2)
library(ggtext)
library(patchwork)

# -- 1. Define your custom order split into two groups of 5
all_pairs <- c(
  "36_35","35_34","34_C","C_A","A_R",
  "A_N1_A_N2","A_N2_A_N3","A_N3_1_N1","1_N1_1_N2","1_N2_1_N3"
)
first5  <- all_pairs[1:5]
second5 <- all_pairs[6:10]

# -- 2. HTML‐styled labels (as before)
label_map <- c(
  "36_35"      = "<span style='color:darkred'>36</span><span style='color:#FFCC99'>_</span><span style='color:darkred'>35</span>",
  "35_34"      = "<span style='color:darkred'>35</span><span style='color:#FFCC99'>_</span><span style='color:darkred'>34</span>",
  "34_C"   = "<span style='color:darkred'>34</span><span style='color:plum'>_</span><span style='color:darkblue'>C</span>",
  "C_A"= "<span style='color:darkblue'>C</span><span style='color:lightblue'>_</span><span style='color:darkblue'>A</span>",
  "A_R" = "<span style='color:darkblue'>A</span><span style='color:lightblue'>_</span><span style='color:darkblue'>R</span>",
  "A_N1_A_N2"  = "<span style='color:darkred'>N1</span><span style='color:#FFCC99'>_</span><span style='color:darkred'>N2</span>",
  "A_N2_A_N3"  = "<span style='color:darkred'>N2</span><span style='color:#FFCC99'>_</span><span style='color:darkred'>N3</span>",
  "A_N3_1_N1"  = "<span style='color:darkred'>N3</span><span style='color:plum'>_</span><span style='color:darkblue'>N1</span>",
  "1_N1_1_N2"  = "<span style='color:darkblue'>N1</span><span style='color:lightblue'>_</span><span style='color:darkblue'>N2</span>",
  "1_N2_1_N3"  = "<span style='color:darkblue'>N2</span><span style='color:lightblue'>_</span><span style='color:darkblue'>N3</span>"
)

# -- 3. Prepare your residuals and Shapiro–Wilk p‐values (unchanged)
library(dplyr)

#Compute residuals for each stack_pair
resid_df <- stack_long_filtered %>%
  group_by(stack_pair) %>%
  do({
    mod <- aov(stack_value ~ car_site * w_position * a_site, data = .)
    tibble(residual = residuals(mod))
  }) %>%
  ungroup()

#Run Shapiro–Wilk on each group
sw_df <- resid_df %>%
  group_by(stack_pair) %>%
  summarise(sw_p = shapiro.test(residual)$p.value) %>%
  ungroup()

# -- 4. Build each of the four panels:

# Panel A: Histograms for first5
hist1 <- ggplot(
  filter(resid_df, stack_pair %in% first5),
  aes(residual)
) +
  geom_histogram(aes(y=..density..), bins=30, fill="grey80", color="white") +
  geom_density(color="steelblue", size=0.6) +
  geom_text(
    data = filter(sw_df, stack_pair %in% first5),
    aes(x=-Inf, y=Inf, label=paste0("SW p=", sw_p)),
    hjust = -0.1, vjust = 1.1, size=2.8, inherit.aes = FALSE
  ) +
  facet_wrap(~stack_pair, nrow=1, ncol=5,
             labeller = labeller(stack_pair = label_map)) +
  theme_minimal(base_size=9) +
  theme(
    strip.text = element_markdown(size=9, face="bold"),
    axis.title = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid.minor = element_blank(),
    panel.spacing = unit(0.3, "lines")
  )

# Panel B: Q–Q for first5
qq1 <- ggplot(
  filter(resid_df, stack_pair %in% first5),
  aes(sample = residual)
) +
  stat_qq(size=0.6) +
  stat_qq_line(color="firebrick", linetype="dashed") +
  facet_wrap(~stack_pair, nrow=1, ncol=5,
             labeller = labeller(stack_pair = label_map)) +
  theme_minimal(base_size=9) +
  theme(
    strip.text = element_markdown(size=9, face="bold"),
    axis.title = element_blank(),
    panel.grid.minor = element_blank(),
    panel.spacing = unit(0.3, "lines")
  )

# Panel C: Histograms for second5
hist2 <- ggplot(
  filter(resid_df, stack_pair %in% second5),
  aes(residual)
) +
  geom_histogram(aes(y=..density..), bins=30, fill="grey80", color="white") +
  geom_density(color="steelblue", size=0.6) +
  geom_text(
    data = filter(sw_df, stack_pair %in% second5),
    aes(x=-Inf, y=Inf, label=paste0("SW p=", sw_p)),
    hjust = -0.1, vjust = 1.1, size=2.8, inherit.aes = FALSE
  ) +
  facet_wrap(~stack_pair, nrow=1, ncol=5,
             labeller = labeller(stack_pair = label_map)) +
  theme_minimal(base_size=9) +
  theme(
    strip.text = element_markdown(size=9, face="bold"),
    axis.title = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid.minor = element_blank(),
    panel.spacing = unit(0.3, "lines")
  )

# Panel D: Q–Q for second5
qq2 <- ggplot(
  filter(resid_df, stack_pair %in% second5),
  aes(sample = residual)
) +
  stat_qq(size=0.6) +
  stat_qq_line(color="firebrick", linetype="dashed") +
  facet_wrap(~stack_pair, nrow=1, ncol=5,
             labeller = labeller(stack_pair = label_map)) +
  theme_minimal(base_size=9) +
  theme(
    strip.text = element_markdown(size=9, face="bold"),
    axis.title = element_blank(),
    panel.grid.minor = element_blank(),
    panel.spacing = unit(0.3, "lines")
  )

# -- 5. Stack all four panels into 4×5
final_fig <- (hist1 / qq1 / hist2 / qq2) +
  plot_annotation(
    title    = "Residual Normality Diagnostics",
    subtitle = "Rows 1 & 3: histograms + density + SW p-values (5 panels each)\nRows 2 & 4: Q–Q plots         (5 panels each)"
  ) &
  theme(
    plot.title    = element_text(size=14, face="bold", hjust=0.5),
    plot.subtitle = element_text(size=9, hjust=0.5),
    plot.margin   = margin(5,5,5,5)
  )

# Render
print(final_fig)


