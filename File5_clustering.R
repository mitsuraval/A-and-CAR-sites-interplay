

##################################################################
# UMAP: Correlation (OR Euclidean Distance)
#################################

library(tidyverse)
library(umap)
library(ggrepel)
library(colorspace)

# ---- STEP 1: Pivot to long format ----
long_data <- combined_hbond_stack_G34 %>%
  pivot_longer(cols = `1`:`30`,
               names_to  = "replicate",
               values_to = "value") %>%
  mutate(
    structure_id = paste(modification, a_site, car_site, w_position, replicate, sep = "_")
  )

# ---- STEP 2: Pivot to wide format ----
umap_input <- long_data %>%
  select(structure_id, bond_name, value) %>%
  pivot_wider(names_from = bond_name, values_from = value) %>%
  select(-y)  # drop problematic column

# ---- STEP 3: Scale the matrix manually ----
umap_matrix <- umap_input %>%
  select(-structure_id) %>%
  as.matrix() %>%
  scale()

# ---- STEP 4a: Compute correlation-distance matrix ----
dist_mat <- 1 - cor(t(umap_matrix), method = "pearson")

# ---- STEP 4b: Prepare a custom config for UMAP ----
custom_config       <- umap.defaults
custom_config$input <- "dist"   # we’re giving it a distance matrix

# ---- STEP 4c: Run UMAP ----
set.seed(42)
# For Euclidean (default), uncomment the next line and comment out the correlation call:
# umap_result <- umap(umap_matrix)
# For correlation distance:
umap_result <- umap(as.matrix(dist_mat), config = custom_config)

umap_coords <- as_tibble(umap_result$layout) %>%
  rename(UMAP1 = V1, UMAP2 = V2)

# ---- STEP 5: Merge metadata ----
umap_data <- bind_cols(
  umap_input %>% select(structure_id),
  umap_coords
) %>%
  separate(structure_id,
           into = c("modification", "a_site", "car_site", "w_position", "replicate"),
           sep  = "_") %>%
  mutate(group_id = paste0(a_site, "-", car_site))

# ---- STEP 6: Re-order facets for a_site ----
# top row: UU, UA, AA; bottom row: CC, UG, AG
umap_data <- umap_data %>%
  mutate(a_site = factor(a_site,
                         levels = c("UU","UA","AA","CC","UG","AG")))

# ---- STEP 7: Color logic ----
a_site_colors <- c(
  "AA" = "red",         "AG" = "blue",       "CC" = "forestgreen",
  "GA" = "purple",      "UA" = "orange",     "UG" = "darkviolet",
  "UU" = "black"
)

umap_data <- umap_data %>%
  mutate(
    group_color = mapply(function(a, car) {
      base <- a_site_colors[[a]]
      if    (car == "GCU") colorspace::darken(base, 0.3)
      else if (car == "CGU") colorspace::lighten(base, 0.3)
      else alpha(base, 0.5)
    }, a_site, car_site)
  )

# ---- STEP 8: Centroids and segment data ----
centroids <- umap_data %>%
  group_by(a_site, car_site, w_position) %>%
  summarise(
    UMAP1 = mean(UMAP1),
    UMAP2 = mean(UMAP2),
    group_color = first(group_color),
    .groups = "drop"
  )

# Compute centroids for GCU vs CGU per (a_site, w_position)
GCU_centroids <- centroids %>%
  filter(car_site == "GCU") %>%
  rename(UMAP1_GCU = UMAP1, UMAP2_GCU = UMAP2)

CGU_centroids <- centroids %>%
  filter(car_site == "CGU") %>%
  rename(UMAP1_CGU = UMAP1, UMAP2_CGU = UMAP2)

segment_data_facet <- inner_join(
  GCU_centroids, CGU_centroids,
  by = c("a_site","w_position")
)

# ---- STEP 9: Final UMAP plot with facets ordered ----
p_umap <- ggplot(umap_data, aes(x = UMAP1, y = UMAP2)) +
  # dashed connectors between GCU & CGU centroids
  geom_segment(
    data = segment_data_facet,
    aes(
      x      = UMAP1_GCU, y      = UMAP2_GCU,
      xend   = UMAP1_CGU, yend   = UMAP2_CGU,
      linetype = w_position
    ),
    color       = "grey50",
    size        = 0.8,
    inherit.aes = FALSE
  ) +
  # scatter points
  geom_point(
    aes(color = car_site, shape = w_position),
    size  = 2.5, alpha = 0.7
  ) +
  # centroids
  geom_point(
    data = centroids,
    aes(x = UMAP1, y = UMAP2, color = car_site, shape = w_position),
    size  = 4, alpha = 1
  ) +
  facet_wrap(~ a_site, ncol = 3) +
  scale_color_manual(
    values = c("GCU" = "darkgreen", "CGU" = "lightgreen"),
    name   = "+1 codon"
  ) +
  scale_shape_manual(
    values = c("wobble" = 16, "wnc" = 17),
    name   = "Decoding geometry",
    labels = c("Wobble (G:U)", "Watson–Crick (G:C)")
  ) +
  scale_linetype_manual(
    values = c("wobble" = "solid", "wnc" = "dashed"),
    guide  = "none"
  ) +
  labs(
    title    = "UMAP of H-bonding & Stacking Fingerprints by A-site",
    subtitle = "Facets = A-site; shapes = decoding geom (● wobble, ▲ WC);\ncolors = +1 codon",
    x = "UMAP 1", y = "UMAP 2"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    strip.text      = element_text(face = "bold", size = 12),
    plot.title      = element_text(face = "bold", size = 18, hjust = 0.5),
    plot.subtitle   = element_text(size = 14, hjust = 0.5),
    legend.position = "bottom"
  )

# ---- Save and view ----
ggsave(
  "Figure_UMAP_all_A_sites.png",
  p_umap,
  width  = 12,
  height = 8,
  dpi    = 300
)
browseURL("Figure_UMAP_all_A_sites.png")








#######################################################################
#UMAP of bonds

# ─── 0. INSTALL & LOAD ────────────────────────────────────────────────────────
if (!requireNamespace("Hmisc", quietly = TRUE)) install.packages("Hmisc")
library(Hmisc)
library(umap)
library(tidyverse)
library(ggrepel)

# ─── BUILD YOUR WIDE MATRIX FOR CORRELATION ───────────────────────────────
wide_mat <- pca_input %>%
  select(-structure_id) %>%   # drop the ID column
  as.matrix()

# ─── GET Pearson CORRELATIONS ─────────────────────────────────────────────
rc    <- rcorr(wide_mat, type = "pearson")
r_mat <- rc$r                   # bond–bond Pearson matrix

# ─── MAKE CORRELATION‐DISTANCE MATRIX ────────────────────────────────────
dist_mat <- 1 - r_mat

# ─── RUN UMAP ON CORRELATION DISTANCES ───────────────────────────────────
config       <- umap.defaults
config$input <- "dist"
set.seed(42)
umap_res     <- umap(as.matrix(dist_mat), config = config)

# ─── Define your bond groups ───────────────────────────────────────────────────
A_site      <- c("a","b","f","g","k","l","m","n","o","p","q")
CAR_site    <- c("d","e","i","j","t","u","v","w","x","z")
Interface   <- c("c","h","r","s")

# ─── Build a named vector of colors ────────────────────────────────────────────
bond_colors <- c(
  # A-site: darkred / #FFCC99 (peach)
  k  = "darkorange2",   a  = "darkorange2",
  l  = "darkorange2",   b  = "darkorange2",
  m  = "darkorange2",   f  = "darkorange2",
  n  = "darkorange2",   g  = "darkorange2",
  o  = "darkorange2",
  p  = "darkorange2",
  q  = "darkorange2",
  
  # CAR-site: darkblue / lightblue
  t  = "deepskyblue3",  d  = "deepskyblue3",
  u  = "deepskyblue3",  e  = "deepskyblue3",
  v  = "deepskyblue3",  i  = "deepskyblue3",
  w  = "deepskyblue3",  j  = "deepskyblue3",
  x  = "deepskyblue3",  z  = "deepskyblue3",
  
  
  # Interface: #74289F (purple) / plum
  r  = "#74289F",   c  = "#74289F",
  s  = "#74289F",   h  = "#74289F"
)

# ─── After you've run the UMAP on dist_mat… ──────────────────────────────────
bond_umap <- tibble(
  bond  = rownames(dist_mat),
  UMAP1 = umap_res$layout[,1],
  UMAP2 = umap_res$layout[,2]
)

# ─── Plot with your custom colors ─────────────────────────────────────────────
ggplot(bond_umap, aes(x = UMAP1, y = UMAP2, color = bond)) +
  geom_point(size = 3, alpha = 0.8) +
  geom_text_repel(aes(label = bond), size = 6) +
  scale_color_manual(values = bond_colors) +
  theme_minimal(base_size = 14) +
  labs(
    title    = "UMAP of Bonds by Correlation Distance",
    subtitle = "A-site (red), CAR-site (blue), Interface (purple)",
    x        = "UMAP 1", y = "UMAP 2"
  )

############


# ── 0) Packages ───────────────────────────────────────────────────────────
library(tidyverse)
library(umap)
library(ggrepel)

# ── 1) Wide matrix (columns = bonds) ──────────────────────────────────────
wide_mat <- pca_input %>%
  select(-structure_id) %>%
  as.matrix()

# Drop constant/empty columns (these cause NA correlations)
const_col <- apply(wide_mat, 2, function(x) sd(x, na.rm = TRUE)) == 0 |
  apply(wide_mat, 2, function(x) all(is.na(x)))
wide_mat  <- wide_mat[, !const_col, drop = FALSE]

# ── 2) Correlation matrix & correlation distance ──────────────────────────
# Choose one:
#   method = "pearson"  (classic)
#   method = "spearman" (rank-robust; often better for noisy/hetero data)
r_mat <- cor(wide_mat, method = "pearson", use = "pairwise.complete.obs")

# Option A (signed):   similarity = r, distance = 1 - r
D <- 1 - r_mat

# (Optional) Option B (sign-agnostic): distance = 1 - |r|
# D <- 1 - abs(r_mat)

# … your code up to building D …

# Clean up numerical issues
diag(D) <- 0
D[!is.finite(D)] <- 0
D <- (D + t(D)) / 2  # force symmetry

# ── 3) UMAP on correlation distances ──────────────────────────────────────
cfg <- umap.defaults
cfg$input <- "dist"   # we're giving UMAP a distance matrix
set.seed(42)

# PASS A MATRIX (not a 'dist' object)
umap_res <- umap(as.matrix(D), config = cfg)


# ── 4) Build a plotting tibble ────────────────────────────────────────────
bond_labels <- colnames(wide_mat)
bond_umap <- tibble(
  bond  = bond_labels,
  UMAP1 = umap_res$layout[, 1],
  UMAP2 = umap_res$layout[, 2]
)

# (Optional) your color map
bond_colors <- c(
  k="darkorange2", a="darkorange2", l="darkorange2", b="darkorange2",
  m="darkorange2", f="darkorange2", n="darkorange2", g="darkorange2",
  o="darkorange2", p="darkorange2", q="darkorange2",
  t="deepskyblue3", d="deepskyblue3", u="deepskyblue3", e="deepskyblue3",
  v="deepskyblue3", i="deepskyblue3", w="deepskyblue3", j="deepskyblue3",
  x="deepskyblue3", z="deepskyblue3",
  r="#74289F", c="#74289F", s="#74289F", h="#74289F"
)

# ── 5) Plot ───────────────────────────────────────────────────────────────
p <- ggplot(bond_umap, aes(UMAP1, UMAP2, color = bond)) +
  geom_point(size = 3, alpha = 0.85) +
  geom_text_repel(aes(label = bond), size = 5, min.segment.length = 0) +
  scale_color_manual(values = bond_colors, na.value = "grey70") +
  theme_minimal(base_size = 14) +
  labs(
    title = "UMAP of Bonds using Correlation Distance",
    subtitle = "Distance = 1 − Spearman correlation (signed). Use 1 − |r| to ignore sign.",
    x = "UMAP 1", y = "UMAP 2", color = "Bond"
  )

print(p)            # show in RStudio
# ggsave("umap_bonds_corrdist.png", p, width = 8, height = 6, dpi = 300)  # optional
