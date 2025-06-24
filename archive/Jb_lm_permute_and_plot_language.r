# script to perform linear model per 100 ms time bin on language-related ERPs with permutation evaluation and FDR correction (Benjamini-Hochberg)
#for only 1 dataset, this is of course all nonsensical!!!

# Required Packages
library(data.table)
library(dplyr)
library(stringr)
library(foreach)
library(doParallel)
library(pbapply)
library(broom)
library(viridis)  # for a nice color scale
library(ggplot2)

script_path <- rstudioapi::getActiveDocumentContext()$path
dir <- dirname(script_path)
setwd(dir)  #set working directory to script path

options(contrasts = c("contr.treatment","contr.Poly"),contr.sum.show.levels = TRUE)

# Function to load a single dataset
load_dataset <- function(file_path) {
  fread(file_path)
}

# Load channel coordinates
define_coords <- function(coords_path) {
  coords <- fread(coords_path, header = FALSE)
  setnames(coords, c("V1", "V2", "V3", "V4"), c("channel", "lat", "sag", "z"))
  coords
}

# Function to load all subjects' data
load_all_subjects <- function(subjects, base_dir, coords) {
  # Initialize parallel backend
  n_cores <- detectCores()
  cl <- makeCluster(n_cores)
  registerDoParallel(cl)
  
  all_subject_data <- foreach(
    s = subjects,
    .packages = c("data.table", "dplyr", "stringr"),
    .export = c("load_dataset")
  ) %dopar% {
    message("Processing subject: ", s)
    # List all files recursively
    all_files <- list.files(base_dir, pattern = "*", recursive = TRUE, full.names = TRUE)
    # Define regex pattern for subject
    pattern <- paste0("^", s, "_[A-Za-z0-9]+_unfold_N\\.csv$")
    # Filter files by basename match
    file_list <- all_files[str_detect(basename(all_files), pattern)]
    # Load and combine
    dfs <- lapply(file_list, load_dataset)
    s_n <- as_tibble(rbindlist(dfs))
    # Process columns
    s_n$channel <- as.character(s_n$channel)
    s_n <- left_join(s_n,coords, by = "channel")
    # Compute laterality
    s_n <- s_n %>%  
      mutate(laterality = ifelse(lat < 0, "left", ifelse(lat > 0, "right", "central")))
    # Compute saggitality
    s_n <- s_n %>%
      mutate(saggitality = case_when(
      sag > 0 & sag <= 0.0714 ~ "frontal",
      sag > 0.0714 ~ "prefrontal",
      sag < 0 & sag >= -0.0929 ~ "posterior",
      sag < -0.0929 ~ "occipital",
      TRUE ~ "central"
      )
    )
   s_n
  }
  
  stopCluster(cl)
  # Combine all subjects
  dat_n <- rbindlist(all_subject_data)
  dat_n
}

# Example usage
subjects <- c('02')#,'03','05','07','09','10','11','12','14','15','16','17','18','19','20','22','23','24','25','26','27','28','29','30','31','33','32')
coords_path <- "r_channel_positions.txt"
coords <- define_coords(coords_path)
base_dir <- paste0(dir,'/eeg/unfold_out/')
dat_n <- load_all_subjects(subjects, base_dir, coords)

dat_n$condition <- as.factor(dat_n$condition)
dat_n$laterality <- as.factor(dat_n$laterality)
dat_n$saggitality <- as.factor(dat_n$saggitality)
# Compute Baseline Factor
baseline_a <- dat_n %>%
  filter(time >= -250, time < 0) %>%
  select(data, subject, condition, laterality, saggitality, mean_target_fixation)

g_baseline <- baseline_a %>%
  group_by(subject, condition, laterality, saggitality, mean_target_fixation) %>%
  summarize(baseline = mean(data, na.rm = TRUE), .groups = "drop")

# Create Time Windows and Join Baseline
create_time_windows <- function(dat, window_size = 100) {
  dat <- dat %>%
    mutate(time_bin = floor(time / window_size) * window_size)
  
  agg <- dat %>%
    group_by(time_bin, laterality, saggitality, condition, mean_target_fixation, subject) %>%
    summarize(data_mean = mean(data, na.rm = TRUE), .groups = "drop")
  
  agg_with_baseline <- inner_join(agg, g_baseline,
                                  by = c("subject", "condition", "laterality", "saggitality", "mean_target_fixation"))
  agg_with_baseline
}


dat_n_windows <- create_time_windows(dat_n) %>% filter(time_bin >= 0 & time_bin <= 1000)



# Model formula
model_formula <- as.formula("data_mean ~ laterality * saggitality * condition * mean_target_fixation * baseline")

# Regression by window
regression_by_window <- function(dat_windows) {
  time_bins <- unique(dat_windows$time_bin)
  results_list <- lapply(time_bins, function(tb) {
    g <- dat_n_windows %>% filter(time_bin == tb)
    model <- lm(model_formula, data = g)
    tidy_out <- broom::tidy(model)
    data.frame(time_bin = tb,
               predictor = tidy_out$term,
               coef = tidy_out$estimate,
               p_value = tidy_out$p.value,
               stringsAsFactors = FALSE)
  })
  do.call(rbind, results_list)
}

regression_results <- regression_by_window(dat_n_windows)
print("Regression Results:")
print(regression_results)

# Block permutation test
block_permutation_test_tstats_fdr <- function(dat_fix_windows,
                                              n_permutations = 2000,
                                              include_baseline = FALSE) {
  library(dplyr)
  
  time_bins <- unique(dat_fix_windows$time_bin)
  
  # 1. Run permutations and get raw p-values
  results_list <- lapply(time_bins, function(tb) {
    g <- dat_fix_windows %>% filter(time_bin == tb)
    
    formula_str <- if (include_baseline) {
      "data_mean ~ laterality * saggitality * condition * mean_target_fixation * baseline"
    } else {
      "data_mean ~ laterality * saggitality * condition * mean_target_fixation"
    }
    formula <- as.formula(formula_str)
    
    # observed t-statistics
    model_obs   <- lm(formula, data = g)
    obs_tstats  <- coef(summary(model_obs))[ , "t value"]
    coef_names  <- names(obs_tstats)
    
    # permutation t-statistics
    perm_tstats <- replicate(n_permutations, {
      g_perm      <- g
      g_perm$data_mean <- sample(g_perm$data_mean)
      model_perm  <- lm(formula, data = g_perm)
      coef(summary(model_perm))[ , "t value"]
    })
    
    # compute empirical (raw) p-values
    p_vals <- sapply(seq_along(obs_tstats), function(k) {
      mean(abs(perm_tstats[k, ]) >= abs(obs_tstats[k]))
    })
    
    data.frame(
      time_bin           = tb,
      predictor          = coef_names,
      observed_tstat     = obs_tstats,
      permutation_p_value= p_vals,
      stringsAsFactors   = FALSE
    )
  })
  
  # 2. Combine and apply FDR within each time_bin
  results_df <- bind_rows(results_list) %>%
    group_by(time_bin) %>%
    mutate(fdr_q_value = p.adjust(permutation_p_value, method = "fdr")) %>%
    ungroup()
  
  return(results_df)
}
# Run block permutation test; for testing, use low n_permutations and include_baseline = FALSE, otherwise it takes quite long
perm_results <- block_permutation_test_tstats_fdr(dat_n_windows, n_permutations = 100, include_baseline = FALSE)

# Write outputs
fwrite(regression_results, "regression_results_nouns.csv")
fwrite(perm_results, "perm_results_nouns.csv")


# Ensure the time_bin column is numeric (if it isn't already)
perm_results <- perm_results %>%
  mutate(time_bin = as.numeric(time_bin))

signif_results = perm_results %>% 
  filter(fdr_q_value <= 0.05 & time_bin >= 0)

find_highest_order <- function(df, p_threshold = 0.05) {
  # Helper: split a predictor string into trimmed parts
  split_predictor <- function(pred) {
    parts <- unlist(strsplit(as.character(pred), "&"))
    trimws(parts)
  }
  
  # Initialize keep column as FALSE
  df$keep <- FALSE
  
  # Process each unique time_bin
  for (tb in unique(df$time_bin)) {
    sub_df <- df[df$time_bin == tb, ]
    
    # Only consider significant rows (p < threshold) and exclude intercept
    sig_df <- sub_df[sub_df$fdr_q_value < p_threshold & sub_df$predictor != "(Intercept)", ]
    if(nrow(sig_df) == 0) next
    
    # Order by number of parts (fewest parts first)
    sig_df$num_parts <- sapply(sig_df$predictor, function(pred) length(split_predictor(pred)))
    sig_df <- sig_df[order(sig_df$num_parts), ]
    
    # For each significant predictor, check if a higher-order predictor exists that includes all its parts.
    for(i in 1:nrow(sig_df)) {
      current_pred <- sig_df$predictor[i]
      cur_parts <- split_predictor(current_pred)
      
      found_superset <- FALSE
      if(i < nrow(sig_df)) {
        for(j in (i+1):nrow(sig_df)) {
          other_pred <- sig_df$predictor[j]
          other_parts <- split_predictor(other_pred)
          # Check if all parts of current are contained in the other AND other has more parts.
          if(length(other_parts) > length(cur_parts) && all(cur_parts %in% other_parts)) {
            found_superset <- TRUE
            break
          }
        }
      }
      sig_df$keep[i] <- !found_superset
    }
    
    # Remove temporary column
    sig_df$num_parts <- NULL
    
    # Update original df for rows in this time_bin.
    for(i in 1:nrow(sig_df)) {
      idx <- which(df$time_bin == tb & df$predictor == sig_df$predictor[i])
      df$keep[idx] <- sig_df$keep[i]
    }
  }
  return(df)
}


# Example usage:
# Assuming perm_results is your permutation results data frame.
# The output will have a new column "keep" (TRUE for highest order, FALSE otherwise).
perm_results_signif <- find_highest_order(perm_results, p_threshold = 0.05)

# Read the permutation results CSV (adjust the filename/path as needed)


perm_results_signif %>% 
  # make sure predictor is a character vector
  mutate(predictor = as.character(predictor)) %>%
  
  # compute the new level ordering:
  { 
    lvl <- unique(.$predictor)
    # count colons in each label
    ncolons <- str_count(lvl, fixed(":"))
    # order by ncolons ascending, then by label alphabetically
    new_levels <- lvl[order(ncolons, lvl)]
    
    # then set predictor as a factor with those levels
    mutate(., predictor = factor(predictor, levels = new_levels))
  } %>%
  mutate(masked_p = ifelse(fdr_q_value > 0.05, NA, fdr_q_value)) %>%  #uncomment if you want plot only significant values (not so god for one dataset)
  mutate(keep = case_when(gsub('baseline','',predictor) != predictor ~ FALSE,
                          TRUE ~ TRUE)) %>% #comment these two lines (253-254) if you want to plot the baseline coefficient, too
  filter(keep == TRUE) %>% #comment this line if you want to plot the basline coefficient, too
  filter(time_bin >= 300 & time_bin <= 700) %>% #time_bins of interest: N4/P3: 300-500, later window 500-800
  drop_na() %>% 
  ggplot(., aes(x = time_bin, y = predictor, fill =fdr_q_value)) +
  geom_tile() +
  scale_fill_viridis(option = "viridis", name = "Permutation p-value (FDR corrected)") 








