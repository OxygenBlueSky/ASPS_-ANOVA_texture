# ASPS Texture Analysis: ANOVA, Post-hoc Tests, and Cohen's d Effect Sizes
#
# Original: 2023.06.06 (Paul)
# Remade: 2025.09.30 (Paul)
# Updated: 2025.10.03 - Added potency decoding functionality
# Updated: 2025.11.03 15:30 CET - Added outer loop for 3 potency filtering scenarios
# Updated: 2025.11.03 15:45 CET - Added Cohen's d effect size calculations
# Updated: 2025.11.26 16:00 CET - Integrated 17-scenario analysis with combined ANOVA + plotting
# Updated: 2025.11.28 15:30 CET - Re-commented with pedagogical style guide
# Updated: 2025.11.28 16:00 CET - Added ratio plot (Stannum/Lactose) with Delta method SE


#===== User settings (edit these) ==========================================

# Set to FALSE to run only the 6-potency scenario (faster for testing/debugging).
# Set to TRUE to run all 17 scenarios (full analysis).
run_all_scenarios <- FALSE


#===== Setup and libraries =================================================

library(openxlsx)
library(readxl)
library(car)
library(here)
library(dplyr)
library(tidyr)
library(ggplot2)
library(plotrix)
library(gridExtra)
library(cowplot)  # For get_legend() in ratio plot.

# Create date string for output filenames (format: YYYYMMDD).
date <- Sys.Date() 
date2 <- gsub("-| |UTC", "", date)

# Create dated output folder; all Excel and PNG files go here.
output_folder <- paste0(date2, "_ASPS_ANOVA_and_plots")
if (!dir.exists(output_folder)) {
  dir.create(output_folder)
  cat(sprintf("Created output folder: %s\n\n", output_folder))
} else {
  cat(sprintf("Using existing output folder: %s\n\n", output_folder))
}


#===== Decoding blinded potencies ==========================================

# Experiments used blinded letter codes (A-F); this lookup table maps them to remedy names.
# Each row = one experiment; columns show which letter was assigned to which remedy.
decoding_file <- "ASPS1-10-decoding table.csv"
decoding_raw <- read.csv(here(decoding_file), skip = 2, stringsAsFactors = FALSE)

colnames(decoding_raw) <- c("experiment_number", "Lactose", "Stannum", "Silicea", 
                            "Sulphur", "Ars.Album", "Mercury")

# Remove the extra header row that came from CSV structure.
decoding_table <- decoding_raw[-1, ]
decoding_table$experiment_number <- as.integer(decoding_table$experiment_number)

cat("============================================================================\n")
cat("DECODING TABLE: Letter codes (A-F) for each remedy by experiment\n")
cat("============================================================================\n")
print(decoding_table)
cat("\n")


# decode_potency: Given experiment number and letter code, returns the remedy name.
# Logic: Look up the experiment row, then search each remedy column to find which one
# contains the matching letter. Returns NA if experiment or letter not found.

decode_potency <- function(experiment_num, potency_code) {
  
  # Get the row for this experiment number.
  exp_row <- decoding_table[decoding_table$experiment_number == experiment_num, ]
  
  if (nrow(exp_row) == 0) {
    return(NA)  # Experiment number not in table.
  }
  
  remedy_cols <- c("Lactose", "Stannum", "Silicea", "Sulphur", "Ars.Album", "Mercury")
  
  # Loop through each remedy column; if the cell value matches the letter code,
  
  # return that column name (= the remedy name).
  for (remedy in remedy_cols) {
    if (exp_row[[remedy]] == potency_code) {
      return(remedy)
    }
  }
  
  return(NA)  # Letter code not found in any column.
}


#===== Load raw data and standardize naming ================================

df <- read.delim(here("ASPS1-10-reduced-data copy.csv"))

# Use only scale=1 observations to avoid pseudoreplication from multiple scales per image.
# The texture analysis software outputs multiple scale levels; we keep only the base level.
df_scale1 <- df[df$scale == 1,]


# AS used series codes (DM, DO, etc.) instead of experiment numbers.
# This lookup maps those codes to experiment numbers 1-5.
series_to_number <- c(
  "DM" = 1,
  "DO" = 2,
  "DQ" = 3,
  "DR" = 4,
  "DS" = 5
)


# parse_substance_group: Standardizes the different naming conventions used by AS and JZ.
# AS format example: "Cress A4 potA" -> need to extract potency letter from end.
# JZ format example: "ASPS 6 A" -> experiment number and potency letter are space-separated.
# Returns a list with: experiment_name (AS or JZ), experiment_number, potency (letter code).

parse_substance_group <- function(substance_group, series_name) {
  
  # Check if it's JZ format (starts with "ASPS").
  if (grepl("^ASPS", substance_group)) {
    # Split on spaces: "ASPS 6 A" -> ["ASPS", "6", "A"]
    parts <- strsplit(trimws(substance_group), " ")[[1]]
    return(list(
      experiment_name = "ASPS_JZ",
      experiment_number = as.integer(parts[2]),  # "6" -> 6
      potency = parts[3]                          # "A"
    ))
  } else {
    # AS format: extract potency letter from end (e.g., "potA" -> "A").
    # regexpr finds position of "pot" followed by letters at string end.
    potency_match <- regmatches(substance_group, 
                                regexpr("pot[A-Z]+$", substance_group))
    
    if (length(potency_match) > 0) {
      potency <- sub("^pot", "", potency_match)  # Remove "pot" prefix.
    } else {
      potency <- "XX"  # Fallback if pattern not found.
    }
    
    return(list(
      experiment_name = "ASPS_AS",
      experiment_number = series_to_number[series_name],  # Look up from series code.
      potency = potency
    ))
  }
}


# Apply parsing to each row, creating standardized columns.
# rowwise() makes mutate() process one row at a time (needed for our parsing function).
df2 <- df_scale1 %>%
  rowwise() %>%
  mutate(
    parsed = list(parse_substance_group(substance_group, series_name)),
    experiment_name = parsed$experiment_name,
    experiment_number = parsed$experiment_number,
    potency = parsed$potency
  ) %>%
  select(-parsed) %>%  # Remove temporary list column.
  ungroup()            # Remove rowwise grouping for normal operations.


# Exclusions

# Exclude Exp 6 chamber 23: dewetting artefact from vaseline caused bad crystallization.
df2 <- df2 %>%
  filter(!(experiment_number == 6 & nr_in_chamber == 23))

cat("Excluded 1 observation due to technical error (Exp 6, chamber 23)\n")
cat("Remaining observations:", nrow(df2), "\n")

# Exclude entire Exp 4: potencies D and E were accidentally mixed during preparation.
df2 <- df2 %>%
  filter(!(experiment_number == 4))

cat("Excluded: Whole Exp 4\n")
cat("Remaining observations:", nrow(df2), "\n")


# Apply decoding to replace letter codes with remedy names

cat("============================================================================\n")
cat("DECODING POTENCIES\n")
cat("============================================================================\n")

# Add new column with actual remedy names by applying decode_potency to each row.
df2 <- df2 %>%
  rowwise() %>%
  mutate(potency_decoded = decode_potency(experiment_number, potency)) %>%
  ungroup()

# Print sample for verification: shows original letter codes alongside decoded names.
cat("Sample of decoded potencies:\n")
print(df2 %>% 
        select(experiment_number, potency, potency_decoded, experiment_name) %>%
        distinct() %>%
        arrange(experiment_number, potency))

cat("\n\nDecoding summary:\n")
cat("Total rows:", nrow(df2), "\n")
cat("Rows with decoded potency:", sum(!is.na(df2$potency_decoded)), "\n")
cat("Rows with missing decoded potency:", sum(is.na(df2$potency_decoded)), "\n")

# Warn if any rows couldn't be decoded (indicates problem with decoding table or data).
if (any(is.na(df2$potency_decoded))) {
  cat("\nWARNING: Some potencies could not be decoded!\n")
  cat("Rows with NA potency_decoded:\n")
  print(df2 %>% 
          filter(is.na(potency_decoded)) %>%
          select(experiment_number, potency, substance_group) %>%
          distinct())
}

cat("\n")

# Scale texture values (×1e6) for numerical stability in ANOVA calculations.
# Very small values can cause precision issues; scaling doesn't affect significance.
df2$maximum_probability <- df2$maximum_probability * 1000000
df2$energy <- df2$energy * 1000000
df2$sum_energy <- df2$sum_energy * 1000000


#===== ANOVA configuration =================================================

# Type III SS requires sum-to-zero contrasts; standard for unbalanced designs
# where group sizes differ. This ensures each effect is tested while controlling
# for all other effects in the model.
options(contrasts = c("contr.sum", "contr.poly"))

# Treat as factors so ANOVA models them as categorical groups, not numeric trends.
df2$potency_decoded <- as.factor(df2$potency_decoded)
df2$experiment_number <- as.factor(df2$experiment_number)


# Define analysis parameters

# Four focal parameters for detailed plotting (pre-selected from preliminary screening).
analysis_vars <- c("cluster_shade", "entropy", "maximum_probability", "kappa")

# All 15 texture parameters for comprehensive Excel ANOVA summary.
all_texture_params <- c(
  "cluster_prominence",
  "cluster_shade",
  "correlation",
  "diagonal_moment",
  "difference_energy",
  "difference_entropy",
  "energy",
  "entropy",
  "inertia",
  "inverse_different_moment",
  "kappa",
  "maximum_probability",
  "sum_energy",
  "sum_entropy",
  "sum_variance"
)

# Hierarchy for ratio calculations: determines which potency is denominator.
# Earlier in list = higher priority to be denominator.
# Lactose (control) preferred, then Stannum (most studied), etc.
potency_hierarchy <- c("Lactose", "Stannum", "Silicea", "Sulphur", "Ars.Album", "Mercury")

# Fixed colors ensure each potency looks the same across all plots for easy comparison.
potency_colors <- c(
  "Lactose" = "#1B9E77",
  "Stannum" = "#D95F02",
  "Silicea" = "#7570B3",
  "Sulphur" = "#E7298A",
  "Ars.Album" = "#66A61E",
  "Mercury" = "#E6AB02"
)


#===== Generate list of 17 potency scenarios ===============================

# Structure: 1 six-potency + 1 three-potency + 15 two-potency combinations = 17 total.
# This allows both full-dataset analysis and focused pairwise comparisons.

cat("\n")
cat("################################################################################\n")
cat("GENERATING ALL 17 POTENCY SCENARIOS\n")
cat("################################################################################\n\n")

all_scenarios <- list()

# Scenario 1: All 6 potencies (full dataset, no filtering).
all_scenarios[[1]] <- list(
  name = "All six potencies",
  suffix = "_6pot",
  potencies = NULL,  # NULL signals "don't filter, use all".
  filename_part = "6pot"
)

# Scenario 2: Three potencies (most commonly used remedies in practice).
all_scenarios[[2]] <- list(
  name = "Stannum, Lactose, Silicea",
  suffix = "_3pot_Stannum-Lactose-Silicea",
  potencies = c("Stannum", "Lactose", "Silicea"),
  filename_part = "3pot_Stannum-Lactose-Silicea"
)

# Scenarios 3-17: All 15 two-potency combinations.
# combn(x, 2) generates all unique pairs from vector x.
all_potencies <- c("Lactose", "Stannum", "Silicea", "Sulphur", "Ars.Album", "Mercury")
two_pot_combos <- combn(all_potencies, 2, simplify = FALSE)  # Returns list of 15 pairs.

for (i in seq_along(two_pot_combos)) {
  combo <- two_pot_combos[[i]]
  all_scenarios[[length(all_scenarios) + 1]] <- list(
    name = paste(combo[1], "vs", combo[2]),
    suffix = paste0("_2pot_", combo[1], "-", combo[2]),
    potencies = combo,
    filename_part = paste0("2pot_", combo[1], "-", combo[2])
  )
}

cat(sprintf("Total scenarios generated: %d\n", length(all_scenarios)))
cat("  1 x 6-potency\n")
cat("  1 x 3-potency\n")
cat("  15 x 2-potency\n\n")


#===== Helper functions ====================================================

# run_anova_summary: Fits Type III ANOVA and extracts p-values.
# Returns dataframe with p-values for: Day (experiment), Potency, and Interaction.
# Falls back to one-way ANOVA (potency only) if data contains only one experiment.

run_anova_summary <- function(data_subset, param_name) {
  
  # Check if parameter column exists in data.
  if (!(param_name %in% colnames(data_subset))) {
    return(data.frame(
      Day = "Error",
      Potency = "Error",
      Interaction = "Error",
      stringsAsFactors = FALSE
    ))
  }
  
  n_experiments <- length(unique(data_subset$experiment_number))
  
  result <- data.frame(
    Day = NA,
    Potency = NA,
    Interaction = NA,
    stringsAsFactors = FALSE
  )
  
  # Wrap in tryCatch to handle cases where model fails (e.g., singular design).
  tryCatch({
    
    if (n_experiments > 1) {
      # Two-way ANOVA with interaction.
      # Tests: (1) potency effect averaged over experiments,
      #        (2) experiment effect averaged over potencies,
      #        (3) whether potency differences vary across experiments.
      formula_str <- paste(param_name, "~ potency_decoded * experiment_number")
      model <- aov(as.formula(formula_str), data = data_subset)
      anova_results <- Anova(model, type = "III")
      
      # Extract p-values from ANOVA table by row name.
      result$Day <- anova_results["experiment_number", "Pr(>F)"]
      result$Potency <- anova_results["potency_decoded", "Pr(>F)"]
      result$Interaction <- anova_results["potency_decoded:experiment_number", "Pr(>F)"]
      
    } else {
      # One-way ANOVA when only one experiment: can only test potency effect.
      formula_str <- paste(param_name, "~ potency_decoded")
      model <- aov(as.formula(formula_str), data = data_subset)
      anova_results <- Anova(model, type = "III")
      
      result$Day <- "NA (single exp)"
      result$Potency <- anova_results["potency_decoded", "Pr(>F)"]
      result$Interaction <- "NA (single exp)"
    }
    
  }, error = function(e) {
    result$Day <- "Error"
    result$Potency <- "Error"
    result$Interaction <- "Error"
  })
  
  return(result)
}


# create_anova_workbook: Generates Excel file with ANOVA results.
# Creates 3 sheets (ALL, AS, JZ) each containing p-values for all 15 texture parameters.
# P-values are color-coded: red < 0.01, orange < 0.05, lilac < 0.10.

create_anova_workbook <- function(df_filtered, scenario) {
  
  wb <- createWorkbook()
  
  # Define cell background colors for significance levels.
  style_red <- createStyle(fgFill = "#FFB3BA")
  style_orange <- createStyle(fgFill = "#FFDFBA")
  style_lilac <- createStyle(fgFill = "#E0BBE4")
  style_header <- createStyle(textDecoration = "bold", fgFill = "#D3D3D3")
  
  # Three data subsets to analyze separately.
  scenarios_excel <- list(
    list(name = "ALL DATA", short_name = "ALL",
         data = df_filtered),
    list(name = "AS ONLY", short_name = "AS",
         data = df_filtered[df_filtered$experiment_name == "ASPS_AS", ]),
    list(name = "JZ ONLY", short_name = "JZ",
         data = df_filtered[df_filtered$experiment_name == "ASPS_JZ", ])
  )
  
  # Loop: For each subset, run ANOVA on all parameters and write results to a sheet.
  for (scen in scenarios_excel) {
    
    addWorksheet(wb, scen$short_name)
    
    # Run ANOVA for each of the 15 texture parameters.
    results_list <- list()
    for (param in all_texture_params) {
      results_list[[param]] <- run_anova_summary(scen$data, param)
    }
    
    # Combine into single dataframe: rows = parameters, cols = Day/Potency/Interaction.
    results_df <- do.call(rbind, results_list)
    results_df <- data.frame(
      Parameter = rownames(results_df),
      results_df,
      stringsAsFactors = FALSE
    )
    
    # Format numeric p-values to 6 decimal places for readability.
    for (col in c("Day", "Potency", "Interaction")) {
      for (row in 1:nrow(results_df)) {
        cell_value <- results_df[row, col]
        if (!is.na(suppressWarnings(as.numeric(cell_value)))) {
          results_df[row, col] <- sprintf("%.6f", as.numeric(cell_value))
        }
      }
    }
    
    # Write title and data table to sheet.
    writeData(wb, scen$short_name,
              paste0("ANOVA Summary: ", scen$name, " - ", scenario$name, " (Exp 4 excluded + outlier removed)"),
              startRow = 1, startCol = 1)
    writeData(wb, scen$short_name, results_df, startRow = 3, rowNames = FALSE)
    
    addStyle(wb, scen$short_name, style_header,
             rows = 3, cols = 1:4, gridExpand = TRUE)
    
    # Apply color-coding to p-value cells based on significance thresholds.
    for (row_idx in 1:nrow(results_df)) {
      for (col_idx in 2:4) {  # Columns 2-4 are Day, Potency, Interaction.
        col_name <- colnames(results_df)[col_idx]
        cell_value <- results_df[row_idx, col_name]
        
        # Skip non-numeric cells (errors, NA messages).
        if (grepl("Error|NA", cell_value)) {
          next
        }
        
        p_val <- suppressWarnings(as.numeric(cell_value))
        
        if (!is.na(p_val)) {
          # Excel row = data row + 3 (accounting for title and header rows).
          excel_row <- row_idx + 3
          excel_col <- col_idx
          
          if (p_val < 0.01) {
            addStyle(wb, scen$short_name, style_red,
                     rows = excel_row, cols = excel_col)
          } else if (p_val < 0.05) {
            addStyle(wb, scen$short_name, style_orange,
                     rows = excel_row, cols = excel_col)
          } else if (p_val < 0.10) {
            addStyle(wb, scen$short_name, style_lilac,
                     rows = excel_row, cols = excel_col)
          }
        }
      }
    }
    
    setColWidths(wb, scen$short_name, cols = 1:4, widths = c(25, 15, 15, 15))
    
    # Add color legend at bottom of sheet.
    legend_row <- nrow(results_df) + 5
    writeData(wb, scen$short_name, "Color Legend:",
              startRow = legend_row, startCol = 1)
    writeData(wb, scen$short_name, "Light lilac = p < 0.10",
              startRow = legend_row + 1, startCol = 1)
    writeData(wb, scen$short_name, "Light orange = p < 0.05",
              startRow = legend_row + 2, startCol = 1)
    writeData(wb, scen$short_name, "Light red = p < 0.01",
              startRow = legend_row + 3, startCol = 1)
    
    addStyle(wb, scen$short_name, style_lilac,
             rows = legend_row + 1, cols = 1)
    addStyle(wb, scen$short_name, style_orange,
             rows = legend_row + 2, cols = 1)
    addStyle(wb, scen$short_name, style_red,
             rows = legend_row + 3, cols = 1)
  }
  
  # Save workbook to output folder.
  output_filename <- file.path(output_folder, 
                               paste0(date2, "_ASPS_ANOVA_", scenario$filename_part, ".xlsx"))
  saveWorkbook(wb, output_filename, overwrite = TRUE)
  
  cat(sprintf("  Excel saved: %s\n", basename(output_filename)))
}


# calculate_ratio: For 2-potency scenarios, computes treatment/control ratio.
# Calculates ratio within each experiment first, then averages across experiments.
# This approach accounts for experiment-to-experiment variability better than
# taking ratio of overall means.

calculate_ratio <- function(df_subset, parameter, potencies) {
  
  # Determine numerator vs denominator using hierarchy.
  # which() returns position in hierarchy vector; lower position = higher priority.
  idx1 <- which(potency_hierarchy == potencies[1])
  idx2 <- which(potency_hierarchy == potencies[2])
  
  # The potency earlier in hierarchy (lower index) becomes denominator (control).
  if (idx1 < idx2) {
    denominator <- potencies[1]
    numerator <- potencies[2]
  } else {
    denominator <- potencies[2]
    numerator <- potencies[1]
  }
  
  if (!(parameter %in% colnames(df_subset))) {
    return(list(numerator = numerator, denominator = denominator, ratio = NA))
  }
  
  experiments <- unique(df_subset$experiment_number)
  experiment_ratios <- numeric(length(experiments))
  
  # Calculate ratio separately for each experiment, then average.
  for (i in seq_along(experiments)) {
    exp <- experiments[i]
    exp_data <- df_subset[df_subset$experiment_number == exp, ]
    
    # Mean of numerator potency in this experiment.
    mean_num_exp <- mean(exp_data[exp_data$potency_decoded == numerator, ][[parameter]], na.rm = TRUE)
    # Mean of denominator potency in this experiment.
    mean_denom_exp <- mean(exp_data[exp_data$potency_decoded == denominator, ][[parameter]], na.rm = TRUE)
    
    experiment_ratios[i] <- mean_num_exp / mean_denom_exp
  }
  
  # Average of per-experiment ratios.
  ratio <- mean(experiment_ratios, na.rm = TRUE)
  
  return(list(
    numerator = numerator,
    denominator = denominator,
    ratio = ratio
  ))
}


# create_plot: Generates 4×2 grid of line plots (4 parameters × 2 experimenters).
# Each plot shows mean ± SE across experiments, with consistent y-axis limits
# for AS vs JZ comparison. For 2-potency scenarios, adds ratio in caption.

create_plot <- function(df_filtered, scenario) {
  
  # Aggregate to mean ± SE per parameter × remedy × experiment × experimenter.
  # std.error() from plotrix calculates SE = SD / sqrt(n).
  summary_data <- df_filtered %>%
    group_by(experiment_name, experiment_number, potency_decoded) %>%
    summarise(
      across(all_of(analysis_vars),
             list(mean = ~mean(.x, na.rm = TRUE),
                  se = ~std.error(.x, na.rm = TRUE)),
             .names = "{.col}_{.fn}"),
      n = n(),
      .groups = "drop"
    )
  
  # Calculate shared y-axis limits per parameter.
  # This ensures AS and JZ plots use same scale for direct visual comparison.
  y_limits <- list()
  for (var in analysis_vars) {
    mean_col <- paste0(var, "_mean")
    se_col <- paste0(var, "_se")
    
    # Find range including error bars, add 5% padding.
    y_min <- min(summary_data[[mean_col]] - summary_data[[se_col]], na.rm = TRUE)
    y_max <- max(summary_data[[mean_col]] + summary_data[[se_col]], na.rm = TRUE)
    
    y_range <- y_max - y_min
    y_limits[[var]] <- c(y_min - 0.05 * y_range, y_max + 0.05 * y_range)
  }
  
  # For 2-potency scenarios, calculate ratios to display in plot captions.
  ratios_as <- list()
  ratios_jz <- list()
  show_ratios <- (!is.null(scenario$potencies) && length(scenario$potencies) == 2)
  
  if (show_ratios) {
    df_as <- df_filtered[df_filtered$experiment_name == "ASPS_AS", ]
    df_jz <- df_filtered[df_filtered$experiment_name == "ASPS_JZ", ]
    
    for (var in analysis_vars) {
      if (nrow(df_as) > 0) {
        ratios_as[[var]] <- calculate_ratio(df_as, var, scenario$potencies)
      }
      if (nrow(df_jz) > 0) {
        ratios_jz[[var]] <- calculate_ratio(df_jz, var, scenario$potencies)
      }
    }
  }
  
  experimenters <- c("ASPS_AS", "ASPS_JZ")
  experimenter_labels <- c("AS (Experiments 1-5)", "JZ (Experiments 6-10)")
  
  plot_list <- list()
  plot_counter <- 1
  
  # Build 8 plots: loop through 4 parameters (rows), then 2 experimenters (columns).
  for (var in analysis_vars) {
    
    mean_col <- paste0(var, "_mean")
    se_col <- paste0(var, "_se")
    
    for (exp_idx in seq_along(experimenters)) {
      
      experimenter <- experimenters[exp_idx]
      experimenter_label <- experimenter_labels[exp_idx]
      
      plot_data <- summary_data %>%
        filter(experiment_name == experimenter)
      
      # Handle case where this experimenter has no data for current scenario.
      if (nrow(plot_data) == 0) {
        plot_list[[plot_counter]] <- ggplot() + 
          theme_void() + 
          ggtitle(paste(experimenter_label, "- No data"))
        plot_counter <- plot_counter + 1
        next
      }
      
      # Dodge width offsets overlapping points horizontally.
      dodge_width <- 0.3
      
      p <- ggplot(plot_data, aes(x = as.numeric(as.character(experiment_number)),
                                 y = .data[[mean_col]],
                                 color = potency_decoded,
                                 group = potency_decoded)) +
        geom_line(linewidth = 0.7, position = position_dodge(width = dodge_width)) +
        geom_errorbar(aes(ymin = .data[[mean_col]] - .data[[se_col]],
                          ymax = .data[[mean_col]] + .data[[se_col]]),
                      width = 0.2, linewidth = 0.5,
                      position = position_dodge(width = dodge_width)) +
        geom_point(size = 2.5, position = position_dodge(width = dodge_width)) +
        scale_x_continuous(breaks = unique(as.numeric(as.character(plot_data$experiment_number)))) +
        scale_y_continuous(limits = y_limits[[var]]) +
        scale_color_manual(values = potency_colors, name = "Remedy") +
        labs(
          title = experimenter_label,
          x = "Experiment Number",
          y = var
        ) +
        theme_bw() +
        theme(
          plot.title = element_text(size = 11, face = "bold"),
          axis.title = element_text(size = 10),
          axis.text = element_text(size = 9),
          legend.title = element_text(size = 9),
          legend.text = element_text(size = 8),
          plot.caption = element_text(size = 8, hjust = 0)
        )
      
      # Add ratio caption for 2-potency scenarios.
      if (show_ratios) {
        ratio_info <- if (experimenter == "ASPS_AS") ratios_as[[var]] else ratios_jz[[var]]
        
        if (!is.null(ratio_info)) {
          caption_text <- sprintf("Ratio: %s/%s = %.3f",
                                  ratio_info$numerator,
                                  ratio_info$denominator,
                                  ratio_info$ratio)
          p <- p + labs(caption = caption_text)
        }
      }
      
      plot_list[[plot_counter]] <- p
      plot_counter <- plot_counter + 1
    }
  }
  
  # Arrange all 8 plots in 4×2 grid.
  grid_plot <- grid.arrange(
    grobs = plot_list,
    nrow = 4,
    ncol = 2,
    top = paste0("ASPS Texture Analysis: ", scenario$name)
  )
  
  # Save to PNG file.
  output_filename <- file.path(output_folder,
                               paste0(date2, "_ASPS_plot_", scenario$filename_part, ".png"))
  
  ggsave(
    filename = output_filename,
    plot = grid_plot,
    width = 30,
    height = 35,
    dpi = 300,
    units = "cm"
  )
  
  cat(sprintf("  Plot saved: %s\n", basename(output_filename)))
}


#===== Main loop: Process scenarios =========================================

# Determine which scenarios to run based on toggle.
if (run_all_scenarios) {
  scenarios_to_run <- seq_along(all_scenarios)
  cat("\n")
  cat("################################################################################\n")
  cat("PROCESSING ALL 17 SCENARIOS\n")
  cat("################################################################################\n\n")
} else {
  scenarios_to_run <- 1  # Only run scenario 1 (all 6 potencies).
  cat("\n")
  cat("################################################################################\n")
  cat("QUICK MODE: Processing only 6-potency scenario (run_all_scenarios = FALSE)\n")
  cat("################################################################################\n\n")
}

for (i in scenarios_to_run) {
  
  scenario <- all_scenarios[[i]]
  
  cat(sprintf("\n========================================================================\n"))
  cat(sprintf("SCENARIO %d/%d: %s\n", i, length(all_scenarios), scenario$name))
  cat(sprintf("========================================================================\n"))
  
  # Filter data to selected potencies (or use all if NULL).
  if (is.null(scenario$potencies)) {
    df_filtered <- df2
    cat("Using all 6 potencies\n")
  } else {
    df_filtered <- df2 %>%
      filter(potency_decoded %in% scenario$potencies)
    cat(sprintf("Filtered to: %s\n", paste(scenario$potencies, collapse = ", ")))
  }
  
  cat(sprintf("Observations: %d\n\n", nrow(df_filtered)))
  
  cat("Creating ANOVA Excel workbook...\n")
  create_anova_workbook(df_filtered, scenario)
  
  cat("Creating plot...\n")
  create_plot(df_filtered, scenario)
  
  cat(sprintf("\nCompleted scenario %d/%d\n", i, length(all_scenarios)))
}

# Report how many files were created.
n_scenarios_run <- length(scenarios_to_run)

cat("\n")
cat("################################################################################\n")
cat("SCENARIO PROCESSING COMPLETED\n")
cat("################################################################################\n")
cat(sprintf("\nScenarios processed: %d of %d\n", n_scenarios_run, length(all_scenarios)))
cat(sprintf("Files created: %d Excel + %d PNG = %d files\n", 
            n_scenarios_run, n_scenarios_run, n_scenarios_run * 2))
cat(sprintf("All outputs saved to: %s/\n", output_folder))


#===== Post-hoc tests: emmeans pairwise comparisons ========================

# Runs on full 6-potency dataset only.
# Computes: (1) pairwise p-values for main potency effect,
#           (2) interaction p-values (does comparison vary across experiments?),
#           (3) Cohen's d effect sizes.

cat("\n\n")
cat("################################################################################\n")
cat("POST HOC TESTS: Full 6-potency dataset only\n")
cat("################################################################################\n\n")

library(emmeans)


# compute_emmeans_contrasts: Fits ANOVA and extracts pairwise comparison results.
# Returns three matrices (all triangular, remedy × remedy):
#   - potency: p-values for main effect pairwise comparisons
#   - interaction: p-values testing whether comparison varies across experiments
#   - cohens_d: effect sizes (mean difference / pooled residual SD)

compute_emmeans_contrasts <- function(data_subset, variable) {
  
  n_experiments <- length(unique(data_subset$experiment_number))
  
  # Fit two-way ANOVA with interaction.
  formula_str <- paste(variable, "~ potency_decoded * experiment_number")
  model <- aov(as.formula(formula_str), data = data_subset)
  
  # Get remedy names from factor levels (determines matrix dimensions).
  remedies <- levels(data_subset$potency_decoded)
  n_remedies <- length(remedies)
  
  # Initialize empty matrices; will fill lower triangle only.
  potency_matrix <- matrix(NA, nrow = n_remedies, ncol = n_remedies,
                           dimnames = list(remedies, remedies))
  interaction_matrix <- matrix(NA, nrow = n_remedies, ncol = n_remedies,
                               dimnames = list(remedies, remedies))
  cohens_d_matrix <- matrix(NA, nrow = n_remedies, ncol = n_remedies,
                            dimnames = list(remedies, remedies))
  
  # Main effect pairwise comparisons: marginal means averaged over experiments.
  # emmeans() calculates estimated marginal means; pairs() does all pairwise contrasts.
  emm_potency <- emmeans(model, ~ potency_decoded)
  pairs_potency <- pairs(emm_potency, adjust = "none")  # No multiple comparison adjustment.
  pairs_summary <- summary(pairs_potency)
  
  # Fill lower triangle of potency matrix with p-values.
  # Loop through each pairwise comparison result.
  for (i in 1:nrow(pairs_summary)) {
    # Parse contrast name (e.g., "Lactose - Stannum") to get remedy names.
    contrast_name <- as.character(pairs_summary$contrast[i])
    parts <- strsplit(contrast_name, " - ")[[1]]
    remedy1 <- trimws(parts[1])
    remedy2 <- trimws(parts[2])
    
    p_val <- pairs_summary$p.value[i]
    
    # Find matrix positions for these remedies.
    idx1 <- which(remedies == remedy1)
    idx2 <- which(remedies == remedy2)
    
    # Place p-value in lower triangle: larger index = row, smaller index = column.
    # This ensures consistent placement regardless of contrast order.
    if (idx1 > idx2) {
      potency_matrix[idx1, idx2] <- p_val
    } else {
      potency_matrix[idx2, idx1] <- p_val
    }
  }
  
  # Cohen's d: standardized effect size = mean difference / residual SD.
  # sigma() extracts residual standard deviation from model (pooled across all groups).
  residual_sd <- sigma(model)
  
  for (i in 1:nrow(pairs_summary)) {
    contrast_name <- as.character(pairs_summary$contrast[i])
    parts <- strsplit(contrast_name, " - ")[[1]]
    remedy1 <- trimws(parts[1])
    remedy2 <- trimws(parts[2])
    
    # estimate = difference in marginal means (remedy1 - remedy2).
    mean_diff <- pairs_summary$estimate[i]
    cohens_d <- mean_diff / residual_sd
    
    # Find matrix positions and place in lower triangle (same logic as p-values).
    idx1 <- which(remedies == remedy1)
    idx2 <- which(remedies == remedy2)
    
    if (idx1 > idx2) {
      cohens_d_matrix[idx1, idx2] <- cohens_d
    } else {
      cohens_d_matrix[idx2, idx1] <- cohens_d
    }
  }
  
  # Interaction test: Does the pairwise comparison vary across experiments?
  # If significant, the main effect may be misleading (effect is inconsistent).
  if (n_experiments > 1) {
    
    tryCatch({
      # Get marginal means separately for each experiment.
      emm_by_exp <- emmeans(model, ~ potency_decoded | experiment_number)
      pairs_by_exp <- pairs(emm_by_exp)  # Pairwise comparisons within each experiment.
      pairs_by_exp_summary <- summary(pairs_by_exp)
      
      # For each remedy pair, test whether the contrast differs across experiments.
      for (i in 1:(n_remedies-1)) {
        for (j in (i+1):n_remedies) {
          remedy1 <- remedies[i]
          remedy2 <- remedies[j]
          
          # Find all rows in pairs_by_exp that match this remedy comparison.
          comparison_name <- paste(remedy1, "-", remedy2)
          matching_rows <- grep(paste0("^", comparison_name, "$"),
                                pairs_by_exp_summary$contrast,
                                fixed = FALSE)
          
          if (length(matching_rows) >= 2) {
            tryCatch({
              # Extract just these contrasts and test if they differ from each other.
              specific_contrasts <- pairs_by_exp[matching_rows]
              pairs_of_contrasts <- pairs(specific_contrasts)  # Contrast of contrasts.
              joint_test <- test(pairs_of_contrasts, joint = TRUE)  # Joint F-test.
              p_val_int <- joint_test$p.value
              
              interaction_matrix[j, i] <- p_val_int
              
            }, error = function(e) {
              interaction_matrix[j, i] <- "Error//pair"
            })
          } else {
            interaction_matrix[j, i] <- "Error/no of exp"
          }
        }
      }
      
    }, error = function(e) {
      interaction_matrix[lower.tri(interaction_matrix)] <- "Error//function"
    })
    
  } else {
    # Can't test interaction with only one experiment.
    interaction_matrix[lower.tri(interaction_matrix)] <- "NA (single exp)"
  }
  
  return(list(
    potency = potency_matrix,
    interaction = interaction_matrix,
    cohens_d = cohens_d_matrix
  ))
}


# Create post-hoc Excel workbook with 6 sheets:
# 3 p-value sheets (ALL, AS, JZ) + 3 Cohen's d sheets (ALL_d, AS_d, JZ_d).

wb_posthoc <- createWorkbook()

# Define cell styles.
style_red <- createStyle(fgFill = "#FFB3BA")
style_orange <- createStyle(fgFill = "#FFDFBA")
style_lilac <- createStyle(fgFill = "#E0BBE4")
style_header <- createStyle(textDecoration = "bold",
                            fgFill = "#D3D3D3",
                            halign = "center",
                            valign = "center")
style_bold <- createStyle(textDecoration = "bold")
style_4dec <- createStyle(numFmt = "0.0000")  # 4 decimal places for p-values.
style_3dec <- createStyle(numFmt = "0.000")   # 3 decimal places for Cohen's d.

# Three data subsets.
scenarios_posthoc <- list(
  list(name = "ALL DATA", short_name = "ALL", data = df2),
  list(name = "AS ONLY", short_name = "AS",
       data = df2[df2$experiment_name == "ASPS_AS", ]),
  list(name = "JZ ONLY", short_name = "JZ",
       data = df2[df2$experiment_name == "ASPS_JZ", ])
)

cat(sprintf("\nLibrary emmeans procedure\n\n"))


# Write p-value sheets (ALL, AS, JZ)
# Each sheet contains triangular matrices for 4 parameters, showing potency and interaction p-values.

for (scenario in scenarios_posthoc) {
  
  cat(sprintf("Processing: %s\n", scenario$name))
  
  addWorksheet(wb_posthoc, scenario$short_name)
  
  current_row <- 1
  
  # Sheet title.
  writeData(wb_posthoc, scenario$short_name,
            paste0("Post Hoc Tests: ", scenario$name, " (Exp 4 excluded + outlier removed)"),
            startRow = current_row, startCol = 1)
  addStyle(wb_posthoc, scenario$short_name, style_bold,
           rows = current_row, cols = 1)
  
  current_row <- current_row + 2
  
  # Loop through 4 focal parameters, building comparison matrix for each.
  for (var in analysis_vars) {
    
    results <- compute_emmeans_contrasts(scenario$data, var)
    
    potency_matrix <- results$potency
    interaction_matrix <- results$interaction
    
    remedies <- rownames(potency_matrix)
    n_remedies <- length(remedies)
    
    # Write header row 1: parameter name + remedy names.
    # Each remedy spans 2 columns (potency p-value + interaction p-value).
    writeData(wb_posthoc, scenario$short_name, var,
              startRow = current_row, startCol = 1)
    
    for (j in 1:n_remedies) {
      writeData(wb_posthoc, scenario$short_name, remedies[j],
                startRow = current_row, startCol = 1 + (j-1)*2 + 1)
      
      # Merge 2 cells for each remedy name.
      mergeCells(wb_posthoc, scenario$short_name,
                 rows = current_row,
                 cols = (1 + (j-1)*2 + 1):(1 + (j-1)*2 + 2))
    }
    
    addStyle(wb_posthoc, scenario$short_name, style_header,
             rows = current_row, cols = 1:(2*n_remedies + 1), gridExpand = TRUE)
    
    current_row <- current_row + 1
    
    # Write header row 2: "potency" and "interaction" labels under each remedy.
    writeData(wb_posthoc, scenario$short_name, "",
              startRow = current_row, startCol = 1)
    
    for (j in 1:n_remedies) {
      writeData(wb_posthoc, scenario$short_name, "potency",
                startRow = current_row, startCol = 1 + (j-1)*2 + 1)
      writeData(wb_posthoc, scenario$short_name, "interaction",
                startRow = current_row, startCol = 1 + (j-1)*2 + 2)
    }
    
    addStyle(wb_posthoc, scenario$short_name, style_header,
             rows = current_row, cols = 1:(2*n_remedies + 1), gridExpand = TRUE)
    
    current_row <- current_row + 1
    
    # Write data rows: Fill lower triangle with p-values.
    # Row i contains comparisons of remedy i vs remedies 1 through i-1.
    # Color-coded by significance: red < 0.01, orange < 0.05, lilac < 0.10.
    for (i in 1:n_remedies) {
      
      # Row label = remedy name.
      writeData(wb_posthoc, scenario$short_name, remedies[i],
                startRow = current_row, startCol = 1)
      
      for (j in 1:n_remedies) {
        
        # Calculate column positions for potency and interaction p-values.
        # Each remedy pair j gets 2 columns: col_pot and col_int.
        col_pot <- 1 + (j-1)*2 + 1
        col_int <- col_pot + 1
        
        # Only fill lower triangle (where column remedy index < row remedy index).
        if (j < i) {
          
          p_pot <- potency_matrix[i, j]
          p_int <- interaction_matrix[i, j]
          
          # Write potency p-value with formatting and significance color.
          if (is.numeric(p_pot) && !is.na(p_pot)) {
            writeData(wb_posthoc, scenario$short_name, p_pot,
                      startRow = current_row, startCol = col_pot)
            
            addStyle(wb_posthoc, scenario$short_name, style_4dec,
                     rows = current_row, cols = col_pot)
            
            # Apply color based on significance threshold.
            if (p_pot < 0.01) {
              addStyle(wb_posthoc, scenario$short_name, style_red,
                       rows = current_row, cols = col_pot, stack = TRUE)
            } else if (p_pot < 0.05) {
              addStyle(wb_posthoc, scenario$short_name, style_orange,
                       rows = current_row, cols = col_pot, stack = TRUE)
            } else if (p_pot < 0.10) {
              addStyle(wb_posthoc, scenario$short_name, style_lilac,
                       rows = current_row, cols = col_pot, stack = TRUE)
            }
          } else if (!is.na(p_pot)) {
            # Write error/NA messages as text.
            writeData(wb_posthoc, scenario$short_name, as.character(p_pot),
                      startRow = current_row, startCol = col_pot)
          }
          
          # Write interaction p-value with same formatting logic.
          if (is.numeric(p_int) && !is.na(p_int)) {
            writeData(wb_posthoc, scenario$short_name, p_int,
                      startRow = current_row, startCol = col_int)
            
            addStyle(wb_posthoc, scenario$short_name, style_4dec,
                     rows = current_row, cols = col_int)
            
            if (p_int < 0.01) {
              addStyle(wb_posthoc, scenario$short_name, style_red,
                       rows = current_row, cols = col_int, stack = TRUE)
            } else if (p_int < 0.05) {
              addStyle(wb_posthoc, scenario$short_name, style_orange,
                       rows = current_row, cols = col_int, stack = TRUE)
            } else if (p_int < 0.10) {
              addStyle(wb_posthoc, scenario$short_name, style_lilac,
                       rows = current_row, cols = col_int, stack = TRUE)
            }
          } else if (!is.na(p_int)) {
            writeData(wb_posthoc, scenario$short_name, as.character(p_int),
                      startRow = current_row, startCol = col_int)
          }
        }
      }
      
      current_row <- current_row + 1
    }
    
    # Add spacing between parameter matrices.
    current_row <- current_row + 2
  }
  
  setColWidths(wb_posthoc, scenario$short_name,
               cols = 1:(2*n_remedies + 1),
               widths = c(15, rep(10, 2*n_remedies)))
  
  # Write legend and interpretation notes.
  legend_row <- current_row + 1
  
  writeData(wb_posthoc, scenario$short_name, "Color Legend:",
            startRow = legend_row, startCol = 1)
  
  writeData(wb_posthoc, scenario$short_name, "Light lilac = p < 0.10",
            startRow = legend_row + 1, startCol = 1)
  addStyle(wb_posthoc, scenario$short_name, style_lilac,
           rows = legend_row + 1, cols = 1)
  
  writeData(wb_posthoc, scenario$short_name, "Light orange = p < 0.05",
            startRow = legend_row + 2, startCol = 1)
  addStyle(wb_posthoc, scenario$short_name, style_orange,
           rows = legend_row + 2, cols = 1)
  
  writeData(wb_posthoc, scenario$short_name, "Light red = p < 0.01",
            startRow = legend_row + 3, startCol = 1)
  addStyle(wb_posthoc, scenario$short_name, style_red,
           rows = legend_row + 3, cols = 1)
  
  writeData(wb_posthoc, scenario$short_name,
            "Note: 'potency' = Main effect p-value (marginal comparison averaged over experiments)",
            startRow = legend_row + 5, startCol = 1)
  
  writeData(wb_posthoc, scenario$short_name,
            "      'interaction' = Does this comparison vary across experiments? (interaction test)",
            startRow = legend_row + 6, startCol = 1)
  
  writeData(wb_posthoc, scenario$short_name,
            "      If interaction p-value is significant, the main effect may be misleading.",
            startRow = legend_row + 7, startCol = 1)
}


# Write Cohen's d sheets (ALL_d, AS_d, JZ_d)
# Each sheet contains triangular matrices of effect sizes for 4 parameters.

cat("\nCreating Cohen's d effect size worksheets\n")

for (scenario in scenarios_posthoc) {
  
  cat(sprintf("Processing Cohen's d: %s\n", scenario$name))
  
  sheet_name_d <- paste0(scenario$short_name, "_d")
  addWorksheet(wb_posthoc, sheet_name_d)
  
  current_row <- 1
  
  # Sheet title.
  writeData(wb_posthoc, sheet_name_d,
            paste0("Cohen's d Effect Sizes: ", scenario$name, " (Exp 4 excluded + outlier removed)"),
            startRow = current_row, startCol = 1)
  addStyle(wb_posthoc, sheet_name_d, style_bold,
           rows = current_row, cols = 1)
  
  current_row <- current_row + 2
  
  # Loop through 4 focal parameters.
  for (var in analysis_vars) {
    
    results <- compute_emmeans_contrasts(scenario$data, var)
    
    cohens_d_matrix <- results$cohens_d
    
    remedies <- rownames(cohens_d_matrix)
    n_remedies <- length(remedies)
    
    # Write header row with parameter name and remedy column labels.
    writeData(wb_posthoc, sheet_name_d, var,
              startRow = current_row, startCol = 1)
    
    for (j in 1:n_remedies) {
      writeData(wb_posthoc, sheet_name_d, remedies[j],
                startRow = current_row, startCol = 1 + j)
    }
    
    addStyle(wb_posthoc, sheet_name_d, style_header,
             rows = current_row, cols = 1:(n_remedies + 1), gridExpand = TRUE)
    
    current_row <- current_row + 1
    
    # Write data rows: Fill lower triangle with Cohen's d values.
    # Positive d = row remedy higher than column remedy.
    for (i in 1:n_remedies) {
      
      writeData(wb_posthoc, sheet_name_d, remedies[i],
                startRow = current_row, startCol = 1)
      
      for (j in 1:n_remedies) {
        
        col_d <- 1 + j
        
        # Only fill lower triangle (j < i).
        if (j < i) {
          
          d_val <- cohens_d_matrix[i, j]
          
          if (is.numeric(d_val) && !is.na(d_val)) {
            writeData(wb_posthoc, sheet_name_d, d_val,
                      startRow = current_row, startCol = col_d)
            
            addStyle(wb_posthoc, sheet_name_d, style_3dec,
                     rows = current_row, cols = col_d)
          } else if (!is.na(d_val)) {
            writeData(wb_posthoc, sheet_name_d, as.character(d_val),
                      startRow = current_row, startCol = col_d)
          }
        }
      }
      
      current_row <- current_row + 1
    }
    
    current_row <- current_row + 2
  }
  
  setColWidths(wb_posthoc, sheet_name_d,
               cols = 1:(n_remedies + 1),
               widths = c(15, rep(10, n_remedies)))
  
  # Write interpretation guide.
  legend_row <- current_row + 1
  
  writeData(wb_posthoc, sheet_name_d, "Cohen's d interpretation:",
            startRow = legend_row, startCol = 1)
  
  writeData(wb_posthoc, sheet_name_d, "Small effect: |d| ≈ 0.2",
            startRow = legend_row + 1, startCol = 1)
  
  writeData(wb_posthoc, sheet_name_d, "Medium effect: |d| ≈ 0.5",
            startRow = legend_row + 2, startCol = 1)
  
  writeData(wb_posthoc, sheet_name_d, "Large effect: |d| ≈ 0.8",
            startRow = legend_row + 3, startCol = 1)
  
  writeData(wb_posthoc, sheet_name_d,
            "Note: Cohen's d calculated using model residual SD (pooled across all groups and experiments)",
            startRow = legend_row + 5, startCol = 1)
  
  writeData(wb_posthoc, sheet_name_d,
            "      Positive d = row remedy has higher value than column remedy",
            startRow = legend_row + 6, startCol = 1)
  
  writeData(wb_posthoc, sheet_name_d,
            "      Negative d = row remedy has lower value than column remedy",
            startRow = legend_row + 7, startCol = 1)
}


# Save post-hoc workbook.

output_filename_posthoc <- file.path(output_folder, 
                                     paste0(date2, "_posthoc_emmeans_cohensd.xlsx"))
saveWorkbook(wb_posthoc, output_filename_posthoc, overwrite = TRUE)

cat(sprintf("\n================================================================================\n"))
cat(sprintf("Post hoc tests exported to: %s\n", basename(output_filename_posthoc)))
cat("================================================================================\n")
cat("File contains 6 sheets:\n")
cat("  P-VALUES: ALL, AS, JZ\n")
cat("  COHEN'S D: ALL_d, AS_d, JZ_d\n")
cat("\nEach sheet has 4 parameters (cluster_shade, entropy, maximum_probability, kappa)\n")
cat("\n================================================================================\n\n")

#===== Ratio plot: Stannum/Lactose across experiments ==========================

# Creates a tall figure with 4 panels (one per focal parameter) showing the
# Stannum/Lactose ratio for each experiment. SE calculated via Delta method.
# NOTE: Delta method assumes large samples for asymptotic normality; with n≈7
# per group, SE estimates are approximate.

cat("################################################################################\n")
cat("RATIO PLOT: Stannum/Lactose across experiments\n")
cat("################################################################################\n\n")


# Filter to Stannum and Lactose only.

df_ratio <- df2 %>%
  filter(potency_decoded %in% c("Stannum", "Lactose"))

cat(sprintf("Filtered to Stannum & Lactose: %d observations\n", nrow(df_ratio)))


# Calculate ratio and SE for each experiment × parameter.
# Delta method for ratio of means (assuming independence):
#   SE_ratio = ratio * sqrt(var_A / (n_A * mean_A^2) + var_B / (n_B * mean_B^2))

ratio_data_list <- list()

for (var in analysis_vars) {
  
  # Get summary stats per experiment × potency.
  summary_by_exp <- df_ratio %>%
    group_by(experiment_number, experiment_name, potency_decoded) %>%
    summarise(
      mean_val = mean(.data[[var]], na.rm = TRUE),
      var_val = var(.data[[var]], na.rm = TRUE),
      n_val = n(),
      .groups = "drop"
    )
  
  # Reshape to wide format: one row per experiment with Stannum and Lactose columns.
  summary_wide <- summary_by_exp %>%
    pivot_wider(
      id_cols = c(experiment_number, experiment_name),
      names_from = potency_decoded,
      values_from = c(mean_val, var_val, n_val)
    )
  
  # Calculate ratio and SE via Delta method.
  summary_wide <- summary_wide %>%
    mutate(
      ratio = mean_val_Stannum / mean_val_Lactose,
      
      # Delta method SE for ratio of independent means.
      # Formula: SE = ratio * sqrt(CV_stan^2 + CV_lac^2)
      # where CV^2 = var / (n * mean^2)
      se_ratio = ratio * sqrt(
        var_val_Stannum / (n_val_Stannum * mean_val_Stannum^2) +
          var_val_Lactose / (n_val_Lactose * mean_val_Lactose^2)
      ),
      
      parameter = var
    )
  
  ratio_data_list[[var]] <- summary_wide
}

# Combine all parameters into single dataframe.
ratio_data <- bind_rows(ratio_data_list)

# Convert experiment_number to numeric for plotting.
ratio_data$experiment_number <- as.numeric(as.character(ratio_data$experiment_number))


# Verification output: Print exp 1 & 6 cluster_shade ratios for manual checking.

cat("\n--- Verification: cluster_shade ratios for Exp 1 and Exp 6 ---\n")
verification_data <- ratio_data %>%
  filter(parameter == "cluster_shade", experiment_number %in% c(1, 6)) %>%
  select(experiment_number, experiment_name, 
         mean_val_Stannum, mean_val_Lactose, 
         var_val_Stannum, var_val_Lactose,
         n_val_Stannum, n_val_Lactose,
         ratio, se_ratio)
print(as.data.frame(verification_data))
cat("----------------------------------------------------------------\n\n")


# X-axis: Keep original experiment numbers (1,2,3,4,5,6,7,8,9,10).
# Exp 4 is excluded from data but kept on axis to show it's missing.
# No position remapping needed.

ratio_data <- ratio_data %>%
  mutate(x_pos = experiment_number)


# Calculate mean ratio per experimenter (for reference lines).
# Also store x-range for each experimenter to limit line extent.

mean_ratios <- ratio_data %>%
  group_by(parameter, experiment_name) %>%
  summarise(
    mean_ratio = mean(ratio, na.rm = TRUE),
    x_min = min(x_pos),
    x_max = max(x_pos),
    .groups = "drop"
  )


# Define colors for AS and JZ.
# Using Sulphur pink for AS, Stannum orange for JZ.
# Lighter versions for reference lines.

color_as <- "#E7298A"
color_jz <- "#D95F02"
color_as_light <- "#F5A9C7"
color_jz_light <- "#FFAA66"


# Set parameter order for consistent panel arrangement.

ratio_data$parameter <- factor(ratio_data$parameter, levels = analysis_vars)
mean_ratios$parameter <- factor(mean_ratios$parameter, levels = analysis_vars)


# Build the plot using individual panels instead of facet_wrap.
# This allows: (1) y-axis label per panel, (2) x-axis on all panels,
# (3) proper segment-based mean lines, (4) section labels with lines.

plot_list <- list()

for (i in seq_along(analysis_vars)) {
  
  var <- analysis_vars[i]
  
  # Subset data for this parameter.
  plot_data <- ratio_data %>% filter(parameter == var)
  mean_data <- mean_ratios %>% filter(parameter == var)
  
  # Get mean ratio values and x-ranges for reference lines.
  mean_as <- mean_data %>% filter(experiment_name == "ASPS_AS")
  mean_jz <- mean_data %>% filter(experiment_name == "ASPS_JZ")
  
  # Build plot for this parameter.
  p <- ggplot(plot_data, aes(x = x_pos, y = ratio, color = experiment_name)) +
    
    # Reference line at ratio = 1 (no effect).
    geom_hline(yintercept = 1, linetype = "dashed", color = "black", linewidth = 0.5) +
    
    # Mean ratio reference lines per experimenter (limited to their x-range).
    geom_segment(
      data = mean_as,
      aes(x = x_min - 0.3, xend = x_max + 0.3, y = mean_ratio, yend = mean_ratio),
      linetype = "dashed", color = color_as_light, linewidth = 0.7,
      inherit.aes = FALSE
    ) +
    geom_segment(
      data = mean_jz,
      aes(x = x_min - 0.3, xend = x_max + 0.3, y = mean_ratio, yend = mean_ratio),
      linetype = "dashed", color = color_jz_light, linewidth = 0.7,
      inherit.aes = FALSE
    ) +
    
    # Points with error bars.
    geom_errorbar(
      aes(ymin = ratio - se_ratio, ymax = ratio + se_ratio),
      width = 0.2, linewidth = 0.5
    ) +
    geom_point(size = 2.5) +
    
    # Color scale.
    scale_color_manual(
      values = c("ASPS_AS" = color_as, "ASPS_JZ" = color_jz),
      labels = c("ASPS_AS" = "AS (Exp 1-5)", "ASPS_JZ" = "JZ (Exp 6-10)"),
      name = "Experimenter"
    ) +
    
    # X-axis: show all experiment numbers including missing exp 4.
    scale_x_continuous(
      breaks = 1:10,
      labels = as.character(1:10),
      limits = c(0.5, 10.5),
      name = "Experiment"
    ) +
    
    # Y-axis: parameter name as label.
    labs(y = var) +
    
    # Theme: clean, modern style.
    theme_bw() +
    theme(
      axis.title.y = element_text(size = 10, face = "bold"),
      axis.title.x = element_text(size = 10),
      axis.text = element_text(size = 9),
      legend.position = "none",  # Will add shared legend later.
      panel.grid.minor = element_blank(),
      plot.margin = margin(t = 5, r = 10, b = 5, l = 10)
    )
  
  # Add section delimiter lines and labels ("AS", "JZ") above data.
  # Calculate y position for labels (top of plot area).
  y_max <- max(plot_data$ratio + plot_data$se_ratio, na.rm = TRUE)
  y_min <- min(plot_data$ratio - plot_data$se_ratio, na.rm = TRUE)
  y_range <- y_max - y_min
  label_y <- y_max + 0.08 * y_range
  line_y <- y_max + 0.02 * y_range
  
  # Expand y-axis to make room for labels.
  p <- p + 
    coord_cartesian(
      ylim = c(y_min - 0.05 * y_range, y_max + 0.15 * y_range),
      clip = "off"
    )
  
  # Add horizontal lines under "AS" and "JZ" labels.
  p <- p +
    annotate("segment", x = 1, xend = 5, y = line_y, yend = line_y,
             color = "gray40", linewidth = 0.5) +
    annotate("segment", x = 6, xend = 10, y = line_y, yend = line_y,
             color = "gray40", linewidth = 0.5) +
    annotate("text", x = 3, y = label_y, label = "AS",
             size = 3.5, fontface = "bold", color = "gray30") +
    annotate("text", x = 8, y = label_y, label = "JZ",
             size = 3.5, fontface = "bold", color = "gray30")
  
  plot_list[[i]] <- p
}


# Create shared legend by extracting from one plot.
# Use guides() to remove the grey background from legend keys.

legend_plot <- ggplot(ratio_data, aes(x = x_pos, y = ratio, color = experiment_name)) +
  geom_point(size = 2.5) +
  scale_color_manual(
    values = c("ASPS_AS" = color_as, "ASPS_JZ" = color_jz),
    labels = c("ASPS_AS" = "AS (Exp 1-5)", "ASPS_JZ" = "JZ (Exp 6-10)"),
    name = "Experimenter"
  ) +
  guides(color = guide_legend(override.aes = list(size = 3))) +
  theme_bw() +
  theme(
    legend.position = "bottom",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9),
    legend.key = element_blank()  # Remove grey background from legend keys.
  )

# Extract legend as grob.
legend_grob <- cowplot::get_legend(legend_plot)

# If cowplot not available, create simple text legend instead.
if (is.null(legend_grob)) {
  # Fallback: add legend to bottom plot.
  plot_list[[4]] <- plot_list[[4]] + theme(legend.position = "bottom")
}


# Combine plots vertically with shared title and legend.

title_grob <- grid::textGrob(
  
  "Treatment Effect Ratio: Stannum vs Lactose\nError bars: SE via Delta method (approximate with n≈7)",
  gp = grid::gpar(fontsize = 12, fontface = "bold"),
  hjust = 0.5
)

# Arrange: title + 4 plots + legend.
combined_plot <- gridExtra::arrangeGrob(
  title_grob,
  plot_list[[1]],
  plot_list[[2]],
  plot_list[[3]],
  plot_list[[4]],
  legend_grob,
  ncol = 1,
  heights = c(0.8, 2.5, 2.5, 2.5, 2.5, 0.5)
)


# Save plot.

output_filename_ratio <- file.path(output_folder,
                                   paste0(date2, "_ASPS_ratio_plot_Stannum_Lactose.png"))

ggsave(
  filename = output_filename_ratio,
  plot = combined_plot,
  width = 14,
  height = 28,
  dpi = 300,
  units = "cm"
)

cat(sprintf("\nRatio plot saved: %s\n", basename(output_filename_ratio)))
cat("================================================================================\n\n")


cat("################################################################################\n")
cat("SCRIPT COMPLETED SUCCESSFULLY\n")
cat("################################################################################\n")
cat(sprintf("\nAll outputs in folder: %s/\n", output_folder))
cat(sprintf("  - %d ANOVA Excel files\n", n_scenarios_run))
cat(sprintf("  - %d scenario plot PNG files\n", n_scenarios_run))
cat(sprintf("  - 1 Post-hoc Excel file\n"))
cat(sprintf("  - 1 Ratio plot PNG file\n"))
cat(sprintf("  TOTAL: %d files\n", n_scenarios_run * 2 + 2))

if (!run_all_scenarios) {
  cat("\nNOTE: Quick mode was used (run_all_scenarios = FALSE).\n")
  cat("Set run_all_scenarios <- TRUE to generate all 17 scenario files.\n")
}