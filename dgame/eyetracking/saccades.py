from typing import Iterable

import numpy as np
import pandas as pd


def compute_saccade_amplitude(start_coords: Iterable, end_coords: Iterable):
    """Compute the amplitude of a saccade from staring and ending coordinates."""
    # Convert to NumPy arrays if necessary
    if not isinstance(start_coords, np.ndarray):
        start_coords = np.array(start_coords, dtype=float)
    if not isinstance(end_coords, np.ndarray):
        end_coords = np.array(end_coords, dtype=float)

    numerator = sum(start_coords * end_coords)
    denominator = np.sqrt(sum(start_coords * start_coords)) * np.sqrt(sum(end_coords * end_coords))
    saccade_amplitude = np.acos(numerator / denominator)
    return saccade_amplitude


def compute_saccade_angles(df: pd.DataFrame) -> pd.DataFrame:
    """Compute radians and degrees of saccades and add to dataframe."""

    # Compute deltas
    df['dx'] = df['norm_pos_x'].diff()
    df['dy'] = df['norm_pos_y'].diff()

    # Compute angles in radians
    df["angles"] = np.arctan2(df["dy"], df["dx"])

    # Convert to degrees for easier interpretation
    df["angles_deg"] = df["angles"] * (360 / np.pi)

    return df
