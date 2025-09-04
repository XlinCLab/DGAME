import argparse
import os

import pandas as pd

from dgame.constants import (CHANNEL_COORDS_FILE, ERP_FIXATION_FILE_SUFFIX,
                             ERP_NOUN_FILE_SUFFIX)
from experiment.load_experiment import Experiment
from experiment.test_subjects import list_subject_files, subject_dirs_dict
from utils.utils import load_csv_list


def main(experiment: str | dict | Experiment) -> Experiment:
    # Initialize DGAME experiment from config
    if not isinstance(experiment, Experiment):
        from dgame.dgame import DGAME
        experiment = DGAME.from_input(experiment)
    
    # Load channel coords file
    channel_coords_file = os.path.join(experiment.input_dir, CHANNEL_COORDS_FILE)  # TODO is this a static file across all subjects? where do we expect it to be found? eeg_dir?
    channel_coords = pd.read_csv(channel_coords_file, names=["channel", "lat", "sag", "z"])
    channel_coords["channel"] = channel_coords["channel"].astype(str)

    # Find per-subject EEG output directories
    per_subject_eeg_outdirs = subject_dirs_dict(
        root_dir = experiment.eeg_outdir,
        subject_regex=experiment.subject_id_regex,
    )
    assert all(len(eeg_outdir) == 1 for eeg_outdir in per_subject_eeg_outdirs.values())
    # Find per-subject unfold_out directories within EEG outdirs
    per_subject_unfold_out_result_dirs = {
        subject_id: os.path.join(subject_eeg_outdirs[0], "unfold_out", "results")
        for subject_id, subject_eeg_outdirs in per_subject_eeg_outdirs.items()
    }
    sorted_subject_ids = sorted(list(per_subject_unfold_out_result_dirs.keys()))

    # Iterate over subjects
    for subject_id in sorted_subject_ids:
        subject_unfold_out_result_dir = per_subject_unfold_out_result_dirs[subject_id]
        erp_noun_files = list_subject_files(subject_unfold_out_result_dir, suffix=ERP_NOUN_FILE_SUFFIX)
        erp_fixation_files = list_subject_files(subject_unfold_out_result_dir, suffix=ERP_FIXATION_FILE_SUFFIX)

        # Load noun data
        erp_noun_data = load_csv_list(erp_noun_files)
        erp_noun_data["channel"] = erp_noun_data["channel"].astype(str)

        # Load fixation data
        erp_fixation_data = load_csv_list(erp_fixation_files)
        erp_fixation_data["channel"] = erp_fixation_data["channel"].astype(str)

        # Merge channel coordinates into both dataframes
        erp_noun_data = erp_noun_data.merge(channel_coords, how="left", on="channel")
        erp_fixation_data = erp_fixation_data.merge(channel_coords, how="left", on="channel")

    return experiment


if __name__ == "__main__":
    parser = argparse.ArgumentParser("Calculate and plot gaze statistics.")
    parser.add_argument('config', help='Path to config.yml file')
    args = parser.parse_args()
    main(os.path.abspath(args.config))
