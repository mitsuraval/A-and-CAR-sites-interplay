

#############################################################

#PLOTTING

#############################################################

# Hbond plots 

# Load required libraries
library(readxl)
library(dplyr)
library(ggplot2)
library(patchwork)

# A-site Plot
A_site_plot <- ggplot(A_filtered_hydrogen_bonds, aes(
  x = a_site,
  y = mean_hydrogen_bonds,
  group = interaction(w_position, car_site),
  color = car_site,
  linetype = w_position
)) +
  geom_point(size = 3) +
  geom_line(linewidth = 1) +
  labs(x = "A-site N1 N2", y = "A site Hydrogen Bonds") +
  theme_minimal() +
  theme(
    axis.text.x  = element_text(angle = 45, hjust = 1, size = 14),
    axis.text.y  = element_text(size = 14),
    axis.title.x = element_text(size = 16, face = "bold"),
    axis.title.y = element_text(size = 16, face = "bold"),
    legend.position = "none"
  ) +
  scale_color_manual(values = c("GCU" = "darkred", "CGU" = "#FFCC99")) +
  scale_linetype_manual(values = c("wobble" = "solid", "wnc" = "dashed"))

# CAR-site Plot
CAR_site_plot <- ggplot(CAR_filtered_hydrogen_bonds, aes(
  x = a_site,
  y = mean_hydrogen_bonds,
  group = interaction(w_position, car_site),
  color = car_site,
  linetype = w_position
)) +
  geom_point(size = 3) +
  geom_line(linewidth = 1) +
  labs(x = "A-site N1 N2", y = "CAR site Hydrogen Bonds") +
  theme_minimal() +
  theme(
    axis.text.x  = element_text(angle = 45, hjust = 1, size = 14),
    axis.text.y  = element_text(size = 14),
    axis.title.x = element_text(size = 16, face = "bold"),
    axis.title.y = element_text(size = 16, face = "bold"),
    legend.position = "none"
  ) +
  scale_color_manual(values = c("GCU" = "darkblue", "CGU" = "lightblue")) +
  scale_linetype_manual(values = c("wobble" = "solid", "wnc" = "dashed"))

# Combine plots side by side
combined_plot <- A_site_plot + CAR_site_plot +
  plot_layout(ncol = 2, guides = "collect")













# Print the combined plot
print(combined_plot)




###########################################################


# STACKING PLOTS

library(ggplot2)
library(dplyr)
library(patchwork)
library(ggtext)

# Define custom order for stack_pair
custom_order <- c(
  "36_35", "35_34", "34_C", "C_A", "A_R",  # Top row
  "A_N1_A_N2", "A_N2_A_N3", "A_N3_1_N1", "1_N1_1_N2", "1_N2_1_N3"       # Bottom row
)

# Define colored facet labels
label_map <- c(
  "36_35" = "<span style='color:darkred'>36</span><span style='color:#FFCC99'>_</span><span style='color:darkred'>35</span>",
  "35_34" = "<span style='color:darkred'>35</span><span style='color:#FFCC99'>_</span><span style='color:darkred'>34</span>",
  "34_C" = "<span style='color:darkred'>34</span><span style='color:plum'>_</span><span style='color:darkblue'>C</span>",
  "C_A" = "<span style='color:darkblue'>C</span><span style='color:lightblue'>_</span><span style='color:darkblue'>A</span>",
  "A_R" = "<span style='color:darkblue'>A</span><span style='color:lightblue'>_</span><span style='color:darkblue'>R</span>",
  "A_N1_A_N2" = "<span style='color:darkred'>N1</span><span style='color:#FFCC99'>_</span><span style='color:darkred'>N2</span>",
  "A_N2_A_N3" = "<span style='color:darkred'>N2</span><span style='color:#FFCC99'>_</span><span style='color:darkred'>N3</span>",
  "A_N3_1_N1" = "<span style='color:darkred'>N3</span><span style='color:plum'>_</span><span style='color:darkblue'>N1</span>",
  "1_N1_1_N2" = "<span style='color:darkblue'>N1</span><span style='color:lightblue'>_</span><span style='color:darkblue'>N2</span>",
  "1_N2_1_N3" = "<span style='color:darkblue'>N2</span><span style='color:lightblue'>_</span><span style='color:darkblue'>N3</span>"
)

# Manually define color for each (stack_pair, car_site) combination
color_map <- data.frame(
  stack_pair = factor(rep(custom_order, each = 2), levels = custom_order),
  car_site = rep(c("GCU", "CGU"), times = length(custom_order)),
  color = c(
    "darkred", "#FFCC99",         # nt36_nt35
    "darkred", "#FFCC99",         # nt35_nt34
    "#74289F", "plum",            # nt34_C1054
    "darkblue", "lightblue",      # C1054_A1196
    "darkblue", "lightblue",      # A1196_R146
    "darkred", "#FFCC99",         # Ant1_Ant2
    "darkred", "#FFCC99",         # Ant2_Ant3
    "#74289F", "plum",            # Ant3_Pnt1
    "darkblue", "lightblue",      # Pnt1_Pnt2
    "darkblue", "lightblue"       # Pnt2_Pnt3
  )
)

# Prepare main data
stack_long_filtered <- stack_long_filtered %>%
  mutate(
    a_site = factor(a_site),
    car_site = factor(car_site),
    w_position = factor(w_position),
    stack_pair = factor(stack_pair, levels = custom_order)
  )

# Summarize mean stacking values
stack_mean <- stack_long_filtered %>%
  group_by(stack_pair, a_site, car_site, w_position) %>%
  summarize(mean_stack = mean(stack_value), .groups = "drop") %>%
  mutate(stack_pair = factor(stack_pair, levels = custom_order)) %>%
  left_join(color_map, by = c("stack_pair", "car_site"))

# Plot
stack_mean_plot <- ggplot(stack_mean, aes(
  x = a_site,
  y = mean_stack,
  group = interaction(w_position, car_site),
  color = color,
  linetype = w_position
)) +
  geom_point(size = 2) +
  geom_line(linewidth = 1) +
  facet_wrap(~stack_pair, ncol = 5, scales = "fixed", labeller = as_labeller(label_map)) +
  labs(
    title = "Mean Stacking Values by A-site, Car-site, and Wobble Position",
    x = "A-site (N1 N2)",
    y = "COG distance (Å)"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    strip.text = element_markdown(face = "bold", size = 13),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  ) +
  scale_color_identity() +
  scale_linetype_manual(values = c("wobble" = "solid", "wnc" = "dashed"))

# Show the plot
print(stack_mean_plot)


############################################################################################

##########violon plots


# Load required libraries
library(ggplot2)
library(dplyr)
library(patchwork)
library(ggtext)

# Define custom order for stack_pair
custom_order <- c(
  "36_35", "35_34", "34_C", "C_A", "A_R",  # Top row
  "A_N1_A_N2", "A_N2_A_N3", "A_N3_1_N1", "1_N1_1_N2", "1_N2_1_N3"       # Bottom row
)

# Manually define color for each (stack_pair, car_site) combination
color_map <- data.frame(
  stack_pair = factor(rep(custom_order, each = 2), levels = custom_order),
  car_site = rep(c("GCU", "CGU"), times = length(custom_order)),
  color = c(
    "darkred", "#FFCC99",         # nt36_nt35
    "darkred", "#FFCC99",         # nt35_nt34
    "#74289F", "plum",            # nt34_C1054
    "darkblue", "lightblue",      # C1054_A1196
    "darkblue", "lightblue",      # A1196_R146
    "darkred", "#FFCC99",         # Ant1_Ant2
    "darkred", "#FFCC99",         # Ant2_Ant3
    "#74289F", "plum",            # Ant3_Pnt1
    "darkblue", "lightblue",      # Pnt1_Pnt2
    "darkblue", "lightblue"       # Pnt2_Pnt3
  )
)

# Define colored facet labels
label_map <- c(
  "36_35" = "<span style='color:darkred'>36</span><span style='color:#FFCC99'>_</span><span style='color:darkred'>35</span>",
  "35_34" = "<span style='color:darkred'>35</span><span style='color:#FFCC99'>_</span><span style='color:darkred'>34</span>",
  "34_C" = "<span style='color:darkred'>34</span><span style='color:plum'>_</span><span style='color:darkblue'>C</span>",
  "C_A" = "<span style='color:darkblue'>C</span><span style='color:lightblue'>_</span><span style='color:darkblue'>A</span>",
  "A_R" = "<span style='color:darkblue'>A</span><span style='color:lightblue'>_</span><span style='color:darkblue'>R</span>",
  "A_N1_A_N2" = "<span style='color:darkred'>N1</span><span style='color:#FFCC99'>_</span><span style='color:darkred'>N2</span>",
  "A_N2_A_N3" = "<span style='color:darkred'>N2</span><span style='color:#FFCC99'>_</span><span style='color:darkred'>N3</span>",
  "A_N3_1_N1" = "<span style='color:darkred'>N3</span><span style='color:plum'>_</span><span style='color:darkblue'>N1</span>",
  "1_N1_1_N2" = "<span style='color:darkblue'>N1</span><span style='color:lightblue'>_</span><span style='color:darkblue'>N2</span>",
  "1_N2_1_N3" = "<span style='color:darkblue'>N2</span><span style='color:lightblue'>_</span><span style='color:darkblue'>N3</span>"
)

# Prepare data and join color information
stack_long_filtered <- stack_long_filtered %>%
  mutate(
    a_site = factor(a_site),
    car_site = factor(car_site),
    w_position = factor(w_position),
    stack_pair = factor(stack_pair, levels = custom_order)
  ) %>%
  left_join(color_map, by = c("stack_pair", "car_site"))

# Generate violin plot
stack_violin_plot <- ggplot(stack_long_filtered, aes(
  x = a_site,
  y = stack_value,
  fill = color
)) +
  geom_violin(scale = "width", alpha = 0.6, color = NA, trim = FALSE) +
  geom_jitter(aes(color = color), size = 0.7, width = 0.2, alpha = 0.5) +
  facet_wrap(
    ~stack_pair,
    ncol = 5,
    scales = "fixed",
    labeller = as_labeller(label_map)
  ) +
  labs(
    title = "Stacking Distributions Across A-site, CAR-site, and Wobble Geometry",
    x = "A-site Codon (N1 N2)",
    y = "COG Distance"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    strip.text = element_markdown(size = 11, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  ) +
  scale_fill_identity() +
  scale_color_identity()

# Show the plot
print(stack_violin_plot)















