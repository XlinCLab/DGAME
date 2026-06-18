# =============================================================================
# extract_instruction_patterns.R
#
# Extracts Director Task instructions from word-based CSVs and derives
# abstract syntactic patterns.
#
# Assumptions about trial segmentation:
#   - The `trial` column is only set for pos %in% c("D","N"); "adjacent"
#     and "word" have trial = 0.
#   - A trial consists of: (a) preamble = adjacent/word rows immediately
#     before the D, (b) D-N core, (c) postamble = adjacent/word rows immediately
#     after the N until the next D or until a large temporal gap.
#   - A "large gap" is defined via GAP_THRESHOLD (default: 2.0 s
#     between consecutive words). This prevents temporally distant words from
#     being incorrectly assigned to a trial.
#
# Output (in OUTPUT_DIR):
#   - instructions_full.csv  : one row per trial with the full word sequence
#   - patterns_abstract.csv  : one row per trial with the slot pattern
#   - pattern_counts.csv     : frequency table of the abstract patterns
#   - pattern_examples.csv   : up to N_EXAMPLES example wordings per pattern
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(purrr)
  library(stringr)
  library(tibble)
})

# ---- Configuration -----------------------------------------------

INPUT_DIR  <- "/media/brilmayer/ieeg_data/experiments/DGame2/preproc/audio"
OUTPUT_DIR <- file.path(INPUT_DIR, "instruction_patterns")
FILE_GLOB  <- "*_words2erp_*_trialtime.csv"   # adjust if needed
GAP_THRESHOLD <- 2.0    # seconds; larger => trial postamble ends
N_EXAMPLES <- 5         # example wordings per pattern

# Direction lemmas that are abstracted to the DIR slot.
DIRECTION_LEMMAS <- c(
  "oben", "unten", "links", "rechts",
  "hoch", "runter", "rauf", "herunter", "hinunter", "hinauf",
  "schraeg", "schräg"
)

# Optional: abstract verbs as VERB slot. Leave empty for the first run.
VERB_LEMMAS <- character(0)
# Example: VERB_LEMMAS <- c("bewege", "schiebe", "raeume", "stelle", "rueck")

dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

# ---- Helper functions ----------------------------------------------------------

read_one <- function(path) {
  df <- suppressMessages(read_csv(
    path, na = c("", "NA"), show_col_types = FALSE
  ))
  df$source_file <- basename(path)
  df$subject_id  <- str_extract(basename(path), "^[0-9]+")
  df
}

#' Assigns each row a trial membership by extending the explicitly
#' annotated D/N trials to their surrounding context.
assign_trial_context <- function(df, gap_threshold = GAP_THRESHOLD) {
  df <- df %>% arrange(time) %>% mutate(row_idx = row_number())
  
  trial_cores <- df %>%
    filter(pos %in% c("D", "N"), !is.na(trial), trial != 0) %>%
    group_by(trial) %>%
    summarise(
      core_start_idx = min(row_idx),
      core_end_idx   = max(row_idx),
      .groups = "drop"
    ) %>%
    arrange(core_start_idx)
  
  df$trial_assigned <- NA_integer_
  
  # Core rows
  for (i in seq_len(nrow(trial_cores))) {
    s <- trial_cores$core_start_idx[i]
    e <- trial_cores$core_end_idx[i]
    df$trial_assigned[s:e] <- trial_cores$trial[i]
  }
  
  # Preamble: backwards from the core
  for (i in seq_len(nrow(trial_cores))) {
    tr <- trial_cores$trial[i]
    s  <- trial_cores$core_start_idx[i]
    j  <- s - 1L
    while (j >= 1L && is.na(df$trial_assigned[j])) {
      gap <- df$time[j + 1L] - df$tmax[j]
      if (is.na(gap) || gap > gap_threshold) break
      df$trial_assigned[j] <- tr
      j <- j - 1L
    }
  }
  
  # Postamble: forwards from the core
  for (i in seq_len(nrow(trial_cores))) {
    tr <- trial_cores$trial[i]
    e  <- trial_cores$core_end_idx[i]
    j  <- e + 1L
    while (j <= nrow(df) && is.na(df$trial_assigned[j])) {
      gap <- df$time[j] - df$tmax[j - 1L]
      if (is.na(gap) || gap > gap_threshold) break
      df$trial_assigned[j] <- tr
      j <- j + 1L
    }
  }
  
  df
}

# Maps a word to its slot type.
classify_slot <- function(text, pos) {
  text_l <- str_to_lower(text)
  is_dir  <- text_l %in% str_to_lower(DIRECTION_LEMMAS)
  is_verb <- length(VERB_LEMMAS) > 0 && text_l %in% str_to_lower(VERB_LEMMAS)
  case_when(
    pos == "D"  ~ "D",
    pos == "N"  ~ "N",
    is_dir      ~ "DIR",
    is_verb     ~ "VERB",
    TRUE        ~ text_l
  )
}

# ---- Main processing ---------------------------------------------------

files <- list.files(INPUT_DIR, pattern = glob2rx(FILE_GLOB),
                    full.names = TRUE, recursive = TRUE)

if (length(files) == 0) {
  stop("Keine Dateien gefunden in ", INPUT_DIR,
       " mit Muster ", FILE_GLOB)
}

message("Verarbeite ", length(files), " Dateien...")

all_trials <- map_dfr(files, function(f) {
  message("  ", basename(f))
  df <- read_one(f) %>% assign_trial_context()
  
  df %>%
    filter(!is.na(trial_assigned)) %>%
    arrange(time) %>%
    group_by(subject_id, source_file, trial_assigned) %>%
    summarise(
      n_words   = n(),
      t_start   = min(time),
      t_end     = max(tmax, na.rm = TRUE),
      duration  = t_end - t_start,
      condition = dplyr::first(stats::na.omit(condition)),
      wording   = paste(text, collapse = " "),
      pattern   = paste(map2_chr(text, pos, classify_slot), collapse = " "),
      has_D     = any(pos == "D"),
      has_N     = any(pos == "N"),
      .groups   = "drop"
    ) %>%
    rename(trial = trial_assigned)
})

# Sanity check
problematic <- all_trials %>% filter(!has_D | !has_N)
if (nrow(problematic) > 0) {
  warning(nrow(problematic), " Trials ohne vollstaendiges D-N-Paar - ",
          "siehe problematic_trials.csv")
  write_csv(problematic, file.path(OUTPUT_DIR, "problematic_trials.csv"))
  all_trials <- all_trials %>% filter(has_D, has_N)
}

# --- Outputs --------------------------------------------------------------

write_csv(
  all_trials %>% select(subject_id, source_file, trial, condition,
                        n_words, duration, wording, pattern),
  file.path(OUTPUT_DIR, "instructions_full.csv")
)

write_csv(
  all_trials %>% select(subject_id, trial, condition, pattern),
  file.path(OUTPUT_DIR, "patterns_abstract.csv")
)

pattern_counts <- all_trials %>%
  count(pattern, sort = TRUE) %>%
  mutate(rel_freq = n / sum(n))
write_csv(pattern_counts, file.path(OUTPUT_DIR, "pattern_counts.csv"))

pattern_examples <- all_trials %>%
  group_by(pattern) %>%
  slice_head(n = N_EXAMPLES) %>%
  ungroup() %>%
  select(pattern, wording, subject_id, trial)
write_csv(pattern_examples, file.path(OUTPUT_DIR, "pattern_examples.csv"))

# ----- Duration statistics --------------------------------------------

# Overall
duration_overall <- all_trials %>%
  summarise(
    n_trials   = n(),
    mean_s     = mean(duration, na.rm = TRUE),
    median_s   = median(duration, na.rm = TRUE),
    sd_s       = sd(duration, na.rm = TRUE),
    min_s      = min(duration, na.rm = TRUE),
    max_s      = max(duration, na.rm = TRUE),
    mean_words = mean(n_words, na.rm = TRUE)
  )

# Per subject
duration_by_subject <- all_trials %>%
  group_by(subject_id) %>%
  summarise(
    n_trials   = n(),
    mean_s     = mean(duration, na.rm = TRUE),
    median_s   = median(duration, na.rm = TRUE),
    sd_s       = sd(duration, na.rm = TRUE),
    min_s      = min(duration, na.rm = TRUE),
    max_s      = max(duration, na.rm = TRUE),
    mean_words = mean(n_words, na.rm = TRUE),
    .groups = "drop"
  )

# Per condition (if present)
duration_by_condition <- all_trials %>%
  filter(!is.na(condition)) %>%
  group_by(condition) %>%
  summarise(
    n_trials   = n(),
    mean_s     = mean(duration, na.rm = TRUE),
    median_s   = median(duration, na.rm = TRUE),
    sd_s       = sd(duration, na.rm = TRUE),
    min_s      = min(duration, na.rm = TRUE),
    max_s      = max(duration, na.rm = TRUE),
    mean_words = mean(n_words, na.rm = TRUE),
    .groups = "drop"
  )

write_csv(duration_overall,     file.path(OUTPUT_DIR, "duration_overall.csv"))
write_csv(duration_by_subject,  file.path(OUTPUT_DIR, "duration_by_subject.csv"))
write_csv(duration_by_condition,file.path(OUTPUT_DIR, "duration_by_condition.csv"))

# ---- Console summary -----------------------------------------

cat("\n========================================\n")
cat("Verarbeitete Trials: ", nrow(all_trials), "\n")
cat("Eindeutige Muster:   ", nrow(pattern_counts), "\n\n")

cat("Instruction duration (s) - gesamt:\n")
print(duration_overall)

cat("\nInstruction duration (s) - pro Subject:\n")
print(duration_by_subject, n = Inf)

if (nrow(duration_by_condition) > 0) {
  cat("\nInstruction duration (s) - pro Condition:\n")
  print(duration_by_condition, n = Inf)
}

cat("\nTop 15 Muster:\n")
print(pattern_counts %>% slice_head(n = 25), n = 25)
cat("\nOutputs in: ", OUTPUT_DIR, "\n")
