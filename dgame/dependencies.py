from dgame.constants import (STEP_DA_KEY, STEP_DB_KEY, STEP_G_KEY, STEP_H_KEY,
                             STEP_I_KEY, STEP_J_KEY)

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
