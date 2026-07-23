# DGAME ANALYSIS STEPS
XDF_EXPORT_DATA_STEP = "xdf.export_audio_and_gaze_times"
WORDS_PREPROCESS_STEP = "words.prepare_words"
WORDS_DESCRIBE_SYNTAX_STEP = "words.describe_syntactic_patterns"
ET_ALIGN_GAZE_AUDIO_STEP = "eyetracking.align_gaze_with_audio"
ET_FIXATION_SACCADES_STEP = "eyetracking.compute_fixation_saccades"
ET_LABEL_TRIALS_STEP = "eyetracking.label_fixation_trials"
ET_GAZE_LANG_STATS_STEP = "eyetracking.gaze_language_timecourse_stats"
ET_PLOT_FIXATION_STEP = "eyetracking.plot_descriptive_fixation_stats"
EEG_PREPROCESS_STEP = "eeg.preprocess_eeg"
EEG_DECONVOLVE_STEP = "eeg.deconvolve_erps"
EEG_RECONSTRUCT_ERPS_STEP = "eeg.reconstruct_erps"
EEG_PLOT_RERPS_STEP = "eeg.plot_rerps"
EEG_REGRESSION_PERMUTATION_STATS_STEP = "eeg.regression_permutation_stats"

# Full pipeline
FULL_DGAME_PIPELINE = [
    # Export data from XDF
    XDF_EXPORT_DATA_STEP,
    WORDS_PREPROCESS_STEP,
    WORDS_DESCRIBE_SYNTAX_STEP,
    # Eyetracking / gaze
    ET_ALIGN_GAZE_AUDIO_STEP,
    ET_FIXATION_SACCADES_STEP,
    ET_LABEL_TRIALS_STEP,
    ET_GAZE_LANG_STATS_STEP,
    ET_PLOT_FIXATION_STEP,
    # EEG
    EEG_PREPROCESS_STEP,
    EEG_DECONVOLVE_STEP,
    EEG_RECONSTRUCT_ERPS_STEP,
    EEG_PLOT_RERPS_STEP,
    EEG_REGRESSION_PERMUTATION_STATS_STEP,
]

# DGAME steps requiring Julia
JULIA_STEPS = {
    EEG_DECONVOLVE_STEP,
    EEG_RECONSTRUCT_ERPS_STEP,
}

# DGAME steps requiring R
R_STEPS = {
    ET_GAZE_LANG_STATS_STEP,
    ET_PLOT_FIXATION_STEP,
    EEG_PLOT_RERPS_STEP,
    EEG_REGRESSION_PERMUTATION_STATS_STEP,
}
