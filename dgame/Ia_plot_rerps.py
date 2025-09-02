import argparse
import os

from dgame.constants import (CHANNEL_FIELD, ERP_FIXATION_FILE_SUFFIX,
                             ERP_NOUN_FILE_SUFFIX)
from experiment.load_experiment import Experiment
from experiment.test_subjects import subject_files_dict
from utils.utils import load_csv_list


def main(experiment: str | dict | Experiment) -> Experiment:
    # Initialize DGAME experiment from config
    if not isinstance(experiment, Experiment):
        from dgame.dgame import DGAME
        experiment = DGAME.from_input(experiment)
    
    # Load channel coords file
    channel_coords = experiment.load_channel_coords()

    # Get selected subject IDs and per-subject files
    subject_erp_noun_files = subject_files_dict(
        dir=os.path.join(experiment.eeg_outdir, "unfold_out", "results"),
        subject_regex=experiment.subject_id_regex,
        suffix=ERP_NOUN_FILE_SUFFIX,
    )
    subject_erp_fixation_files = subject_files_dict(
        dir=os.path.join(experiment.eeg_outdir, "unfold_out", "results"),
        subject_regex=experiment.subject_id_regex,
        suffix=ERP_FIXATION_FILE_SUFFIX,
    )
    subject_ids = noun_subject_ids = sorted(list(subject_erp_noun_files.keys()))
    fixation_subject_ids = sorted(list(subject_erp_fixation_files.keys()))
    # Make sure the same subject IDs are found for both types
    assert fixation_subject_ids == noun_subject_ids

    # Iterate over subjects
    for subject_id in subject_ids:
        erp_noun_files = subject_erp_noun_files[subject_id]
        erp_fixation_files = subject_erp_fixation_files[subject_id]

        # Load noun data
        erp_noun_data = load_csv_list(erp_noun_files)
        erp_noun_data[CHANNEL_FIELD] = erp_noun_data[CHANNEL_FIELD].astype(str)

        # Load fixation data
        erp_fixation_data = load_csv_list(erp_fixation_files)
        erp_fixation_data[CHANNEL_FIELD] = erp_fixation_data[CHANNEL_FIELD].astype(str)

        # Merge channel coordinates into both dataframes
        erp_noun_data = erp_noun_data.merge(channel_coords, how="left", on=CHANNEL_FIELD)
        erp_fixation_data = erp_fixation_data.merge(channel_coords, how="left", on=CHANNEL_FIELD)

    return experiment


if __name__ == "__main__":
    parser = argparse.ArgumentParser("Calculate and plot gaze statistics.")
    parser.add_argument('config', help='Path to config.yml file')
    args = parser.parse_args()
    main(os.path.abspath(args.config))
