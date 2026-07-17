# DGAME ANALYSIS STEPS
STEP_A_KEY = "A_export_audio_and_et_times"
STEP_B_KEY = "words.prepare_words"
STEP_CA_KEY = "eyetracking.align_gaze_with_audio"
STEP_CB_KEY = "eyetracking.compute_fixation_saccades"
STEP_CC_KEY = "eyetracking.label_fixation_trials"
STEP_DA1_KEY = "Da_compute_trialtime"
STEP_DA2_KEY = "Da_gaze_language_timecourse_stats"
STEP_DB_KEY = "eyetracking.plot_descriptive_fixation_stats"
STEP_E_KEY = "words.describe_syntactic_patterns"
STEP_F_KEY = "F_preproc_EEG"
STEP_G_KEY = "G_deconvolution_ERPs"
STEP_H_KEY = "H_reconstruct_ERPs"
STEP_I_KEY = "I_plot_rERPs"
STEP_J_KEY = "J_lm_permute_and_plot_fixations_and_language"

# Full pipeline
FULL_DGAME_PIPELINE = [
    STEP_A_KEY,
    STEP_B_KEY,
    STEP_CA_KEY,
    STEP_CB_KEY,
    STEP_CC_KEY,
    STEP_DA1_KEY,
    STEP_DA2_KEY,
    STEP_DB_KEY,
    STEP_E_KEY,
    STEP_F_KEY,
    STEP_G_KEY,
    STEP_H_KEY,
    STEP_I_KEY,
    STEP_J_KEY,
]

# DGAME steps requiring Julia
JULIA_STEPS = {
    STEP_G_KEY,
    STEP_H_KEY,
}

# DGAME steps requiring R
R_STEPS = {
    STEP_DA2_KEY,
    STEP_DB_KEY,
    STEP_I_KEY,
    STEP_J_KEY,
}
