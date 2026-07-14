# DGAME ANALYSIS STEPS
STEP_A_KEY = "A_export_audio_and_et_times"
STEP_B_KEY = "B_prepare_words"
STEP_CA_KEY = "Ca_preproc_et_data"
STEP_CB_KEY = "Cb_preproc_fixations"
STEP_CC_KEY = "Cc_prepare_fixations"
STEP_DA_KEY = "Da_gaze_stats"
STEP_DB_KEY = "Db_plot_descriptive_fixation"
STEP_E_KEY = "E_describe_syntactic_patterns_from_audio_instructions"
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
    STEP_DA_KEY,
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
    STEP_DA_KEY,
    STEP_DB_KEY,
    STEP_I_KEY,
    STEP_J_KEY,
}
