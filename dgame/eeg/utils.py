import numpy as np
import pandas as pd

from dgame.eeg import (LATERAL_INPUT_FIELD, LATERALITY_FIELD,
                       SAGGITAL_INPUT_FIELD, SAGGITALITY_FIELD)


def annotate_laterality_and_saggitality(df: pd.DataFrame) -> pd.DataFrame:
    """Annotate pandas DataFrame with laterality and saggitality labels."""
    if SAGGITAL_INPUT_FIELD not in df.columns:
        raise ValueError(f"DataFrame missing '{SAGGITAL_INPUT_FIELD}' column")
    if LATERAL_INPUT_FIELD not in df.columns:
        raise ValueError(f"DataFrame missing '{LATERAL_INPUT_FIELD}' column")

    # Compute laterality
    df[LATERALITY_FIELD] = np.where(
        df[LATERAL_INPUT_FIELD] < 0, "left",
        np.where(df[LATERAL_INPUT_FIELD] > 0, "right", "central")
    )

    # Compute saggitality
    sag_conditions = [
        (df[SAGGITAL_INPUT_FIELD] > 0) & (df[SAGGITAL_INPUT_FIELD] <= 0.0714),      # frontal
        (df[SAGGITAL_INPUT_FIELD] > 0.0714),                                        # prefrontal
        (df[SAGGITAL_INPUT_FIELD] < 0) & (df[SAGGITAL_INPUT_FIELD] >= -0.0929),     # posterior
        (df[SAGGITAL_INPUT_FIELD] < -0.0929)                                        # occipital
        # elsewhere condition                                                       # central
    ]
    sag_labels = ["frontal", "prefrontal", "posterior", "occipital"]
    df[SAGGITALITY_FIELD] = np.select(sag_conditions, sag_labels, default="central")

    return df
