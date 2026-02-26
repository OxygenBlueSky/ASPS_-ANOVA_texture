# Correlation Matrix Plots: HOM2 + ASPS
# Produces correlation matrices for both datasets in a dated output folder:
#   - HOM2: 16×16 (15 texture + EVAP), split by ALL / LAB_1 / LAB_2
#   - ASPS texture: 16×16 (15 texture + evaporation), split by ALL / AS / JZ
#   - ASPS combined: 19×19 (15 texture + 3 fractal + evaporation), split by ALL / AS / JZ


#===== Setup and libraries =================================================

library(readxl)
library(dplyr)
library(ggplot2)
library(GGally)
library(here)


#===== Dated output folder =================================================

date2 <- format(Sys.Date(), "%Y%m%d")
output_folder <- paste0(date2, "_correlation_plots")

if (!dir.exists(output_folder)) {
  dir.create(output_folder)
  cat(sprintf("Created output folder: %s\n\n", output_folder))
} else {
  cat(sprintf("Using existing output folder: %s\n\n", output_folder))
}


#===== Shared panel functions ==============================================

# upper_cor: Upper triangle - Pearson r with green color gradient
upper_cor <- function(data, mapping, ...) {
  x <- eval_data_col(data, mapping$x)
  y <- eval_data_col(data, mapping$y)
  cor_val <- cor(x, y, use = "complete.obs")

  abs_cor <- abs(cor_val)
  red_component <- round(144 * (1 - abs_cor))
  green_component <- round(238 - 138 * abs_cor)
  blue_component <- round(144 * (1 - abs_cor))
  color_hex <- sprintf("#%02X%02X%02X", red_component, green_component, blue_component)

  ggplot(data, mapping) +
    annotate("text", x = 0.5, y = 0.5,
             label = sprintf("%.2f", cor_val),
             color = color_hex, size = 6) +
    xlim(0, 1) + ylim(0, 1) +
    theme_void()
}

# lower_scatter: Lower triangle - scatterplots
lower_scatter <- function(data, mapping, ...) {
  ggplot(data, mapping) +
    geom_point(alpha = 0.2, size = 0.5, color = "gray40") +
    theme_minimal() +
    theme(
      panel.grid = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank()
    )
}

# diag_label: Diagonal - parameter names
diag_label <- function(data, mapping, ...) {
  var_name <- rlang::as_label(mapping$x)
  display_name <- gsub("_", "\n", var_name)

  ggplot(data, mapping) +
    annotate("text", x = 0.5, y = 0.5,
             label = display_name,
             size = 2.5, fontface = "bold", lineheight = 0.8) +
    xlim(0, 1) + ylim(0, 1) +
    theme_void()
}


#===== Shared plotting function ============================================

make_corr_plot <- function(data, param_cols, title, output_path, size_cm = 40) {

  cat(sprintf("  Generating: %s (%d observations)...\n", title, nrow(data)))

  df_ordered <- data[, param_cols]

  plot_corr <- ggpairs(
    df_ordered,
    upper = list(continuous = upper_cor),
    lower = list(continuous = lower_scatter),
    diag = list(continuous = diag_label)
  ) +
    theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      strip.text = element_blank(),
      panel.spacing = unit(0.1, "lines")
    ) +
    labs(title = title)

  ggsave(
    filename = output_path,
    plot = plot_corr,
    width = size_cm, height = size_cm, dpi = 300, units = "cm"
  )

  cat(sprintf("  Saved: %s\n", basename(output_path)))
}


# parse_fractal_filename: Extract series and chamber from fractal filename
# JZ format: "H120220909GCA-01-20220924..." -> series=GCA, chamber=01
# AS format: "P20160416DM-01-20160824..."   -> series=DM,  chamber=01
# (Reused from main ASPS script)

parse_fractal_filename <- function(filename) {

  # Try JZ format first (GC[A-E]-[0-9]+)
  if (grepl("GC[A-E]", filename)) {
    match <- regmatches(filename, regexpr("GC[A-E]-[0-9]+", filename))
    parts <- strsplit(match, "-")[[1]]
    return(list(
      experiment_name = "JZ",
      series_name = parts[1],
      nr_in_chamber = as.integer(parts[2])
    ))
  }

  # Try AS format (P[0-9]+D[MOPQRS]-[0-9]+)
  if (grepl("P[0-9]+D[MOPQRS]", filename)) {
    match <- regmatches(filename, regexpr("P[0-9]+D[MOPQRS]-[0-9]+", filename))
    clean <- sub("P[0-9]*", "", match)
    parts <- strsplit(clean, "-")[[1]]
    return(list(
      experiment_name = "AS",
      series_name = parts[1],
      nr_in_chamber = as.integer(parts[2])
    ))
  }

  return(list(experiment_name = NA, series_name = NA, nr_in_chamber = NA))
}


#===== HOM2 DATA ===========================================================

cat("============================================================================\n")
cat("HOM2 DATA\n")
cat("============================================================================\n\n")

df_hom2_raw <- read_xls(here("HOM2data.xls"))
cat("Loaded HOM2 data:", nrow(df_hom2_raw), "rows\n")

# Filter to scale=1 (avoid pseudoreplication from multiple scales per image)
df_hom2 <- df_hom2_raw[df_hom2_raw$scale == 1, ]
cat("After scale=1 filter:", nrow(df_hom2), "rows\n")

# Keep only rows where MEX_VER2 is not blank and is between 1-20
df_hom2 <- df_hom2 %>%
  filter(!is.na(MEX_VER2), MEX_VER2 >= 1, MEX_VER2 <= 20)
cat("After MEX_VER2 filter (1-20, non-blank):", nrow(df_hom2), "rows\n")

# Exclude Autoclave = 1 and blanks
df_hom2 <- df_hom2 %>%
  filter(!is.na(Autoclave), Autoclave != 1)
cat("After excluding Autoclave=1 and blanks:", nrow(df_hom2), "rows\n")

# Exclude CellPhone = 1 and blanks
df_hom2 <- df_hom2 %>%
  filter(!is.na(CellPhone), CellPhone != 1)
cat("After excluding CellPhone=1 and blanks:", nrow(df_hom2), "rows\n")

cat(sprintf("\nHOM2 final dataset: %d observations\n\n", nrow(df_hom2)))

# HOM2 parameter order: 15 texture + EVAP
hom2_params <- c(
  "cluster_shade", "diagonal_moment", "maximum_probability", "sum_energy",
  "cluster_prominence", "entropy", "kappa", "energy", "correlation",
  "difference_energy", "difference_entropy", "inertia",
  "inverse_different_moment", "sum_entropy", "sum_variance",
  "EVAP"
)

# HOM2 subsets: ALL, LAB=1, LAB=2
hom2_subsets <- list(
  list(label = "ALL",   data = df_hom2),
  list(label = "LAB_1", data = df_hom2 %>% filter(LAB == 1)),
  list(label = "LAB_2", data = df_hom2 %>% filter(LAB == 2))
)

for (s in hom2_subsets) {
  make_corr_plot(
    data = s$data,
    param_cols = hom2_params,
    title = paste("Correlation Matrix: HOM2", s$label, "(Texture Parameters)"),
    output_path = file.path(output_folder,
                            paste0("HOM2_texture_correlation_", s$label, ".png"))
  )
}

cat("\nHOM2 plots done.\n\n")


#===== ASPS DATA ===========================================================

cat("============================================================================\n")
cat("ASPS DATA\n")
cat("============================================================================\n\n")

df_asps_raw <- read.delim(here("ASPS1-10-reduced-data copy.csv"))
cat("Loaded ASPS data:", nrow(df_asps_raw), "rows\n")

# Filter to scale=1
df_asps <- df_asps_raw[df_asps_raw$scale == 1, ]
cat("After scale=1 filter:", nrow(df_asps), "rows\n")

# Exclude Exp 4 (series DR): potencies D and E were accidentally mixed
initial_n <- nrow(df_asps)
df_asps <- df_asps %>%
  filter(series_name != "DR")
cat(sprintf("Excluded Exp 4 / DR (mixed potencies): %d rows removed\n",
            initial_n - nrow(df_asps)))

# Exclude Exp 6 chamber 23 (series GCA): dewetting artefact
initial_n <- nrow(df_asps)
df_asps <- df_asps %>%
  filter(!(series_name == "GCA" & nr_in_chamber == 23))
cat(sprintf("Excluded Exp 6 / GCA chamber 23 (technical error): %d rows removed\n",
            initial_n - nrow(df_asps)))

# Exclude reference plates (position swap between Exp 1 and Exp 2-10)
# Exp 1 (DM): reference at position 4, so exclude position 8
# Exp 2-10: reference at position 8, so exclude position 4
initial_n <- nrow(df_asps)
df_asps <- df_asps %>%
  filter(!(series_name == "DM" & nr_in_chamber == 8)) %>%
  filter(!(series_name != "DM" & nr_in_chamber == 4))
cat(sprintf("Excluded reference plates: %d rows removed\n",
            initial_n - nrow(df_asps)))

# Derive experiment_name from series_name
df_asps <- df_asps %>%
  mutate(experiment_name = ifelse(
    series_name %in% c("DM", "DO", "DQ", "DR", "DS"), "AS", "JZ"
  ))

cat(sprintf("\nASPS final dataset: %d observations\n", nrow(df_asps)))
cat(sprintf("  AS: %d,  JZ: %d\n\n",
            sum(df_asps$experiment_name == "AS"),
            sum(df_asps$experiment_name == "JZ")))

# ASPS parameter order: 15 texture + evaporation_duration
asps_params <- c(
  "cluster_shade", "diagonal_moment", "maximum_probability", "sum_energy",
  "cluster_prominence", "entropy", "kappa", "energy", "correlation",
  "difference_energy", "difference_entropy", "inertia",
  "inverse_different_moment", "sum_entropy", "sum_variance",
  "evaporation_duration"
)

# ASPS subsets: ALL, AS, JZ
asps_subsets <- list(
  list(label = "ALL", data = df_asps),
  list(label = "AS",  data = df_asps %>% filter(experiment_name == "AS")),
  list(label = "JZ",  data = df_asps %>% filter(experiment_name == "JZ"))
)

for (s in asps_subsets) {
  make_corr_plot(
    data = s$data,
    param_cols = asps_params,
    title = paste("Correlation Matrix: ASPS", s$label, "(Texture Parameters)"),
    output_path = file.path(output_folder,
                            paste0("ASPS_texture_correlation_", s$label, ".png"))
  )
}

cat("\nASPS texture-only plots done.\n\n")


#===== ASPS COMBINED (TEXTURE + FRACTAL) ===================================

cat("============================================================================\n")
cat("ASPS COMBINED: TEXTURE + FRACTAL\n")
cat("============================================================================\n\n")

# Load fractal data
df_frac_raw <- read.delim(here("20251217_fractal_data_ASPS_1-10 copy.csv"))
cat("Loaded fractal data:", nrow(df_frac_raw), "rows\n")

# Rename cryptic FracLac column names to readable parameter names
colnames(df_frac_raw)[6] <- "db_mean"
colnames(df_frac_raw)[17] <- "dm_mean"
colnames(df_frac_raw)[26] <- "dx_mean"

# Parse filenames to extract series_name and nr_in_chamber
df_frac <- data.frame(
  filename = df_frac_raw[, 1],
  db_mean = df_frac_raw$db_mean,
  dm_mean = df_frac_raw$dm_mean,
  dx_mean = df_frac_raw$dx_mean,
  stringsAsFactors = FALSE
)

df_frac <- df_frac %>%
  rowwise() %>%
  mutate(
    parsed = list(parse_fractal_filename(filename)),
    series_name = parsed$series_name,
    nr_in_chamber = parsed$nr_in_chamber
  ) %>%
  select(-parsed, -filename) %>%
  ungroup()

cat(sprintf("Parsed fractal data: %d rows with valid series/chamber\n",
            sum(!is.na(df_frac$series_name))))

# Inner join with existing df_asps (already filtered, exclusions applied)
df_combined <- df_asps %>%
  inner_join(df_frac, by = c("series_name", "nr_in_chamber"))

cat(sprintf("\nAfter inner join: %d matched rows\n", nrow(df_combined)))
cat(sprintf("  AS: %d,  JZ: %d\n\n",
            sum(df_combined$experiment_name == "AS"),
            sum(df_combined$experiment_name == "JZ")))

# Combined parameter order: 15 texture + 3 fractal + evaporation_duration (19 total)
combined_params <- c(
  "cluster_shade", "diagonal_moment", "maximum_probability", "sum_energy",
  "cluster_prominence", "entropy", "kappa", "energy", "correlation",
  "difference_energy", "difference_entropy", "inertia",
  "inverse_different_moment", "sum_entropy", "sum_variance",
  "db_mean", "dm_mean", "dx_mean",
  "evaporation_duration"
)

# Combined subsets: ALL, AS, JZ
combined_subsets <- list(
  list(label = "ALL", data = df_combined),
  list(label = "AS",  data = df_combined %>% filter(experiment_name == "AS")),
  list(label = "JZ",  data = df_combined %>% filter(experiment_name == "JZ"))
)

for (s in combined_subsets) {
  make_corr_plot(
    data = s$data,
    param_cols = combined_params,
    title = paste("Correlation Matrix: ASPS", s$label, "(Texture + Fractal)"),
    output_path = file.path(output_folder,
                            paste0("ASPS_texture_fractal_correlation_", s$label, ".png")),
    size_cm = 50
  )
}

cat("\nASPS combined (texture + fractal) plots done.\n")


#===== Summary =============================================================

cat("\n============================================================================\n")
cat("ALL DONE\n")
cat("============================================================================\n")
cat(sprintf("9 correlation matrix PNGs saved to: %s/\n", output_folder))
