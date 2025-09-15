import argparse
import os

from dgame.constants import (CHANNEL_FIELD, ERP_FIXATION_FILE_SUFFIX,
                             ERP_NOUN_FILE_SUFFIX)
from dgame.Ja_lm_permute_and_plot_fixations import \
    annotate_laterality_and_saggitality
from experiment.load_experiment import Experiment
from experiment.test_subjects import list_subject_files
from utils.utils import load_csv_list


def main(experiment: str | dict | Experiment) -> Experiment:
    # Initialize DGAME experiment from config
    if not isinstance(experiment, Experiment):
        from dgame.dgame import DGAME
        experiment = DGAME.from_input(experiment)

    # Find per-subject EEG output directories
    per_subject_eeg_outdirs = experiment.get_subject_dirs_dict(experiment.eeg_outdir)
    assert all(len(eeg_outdir) == 1 for eeg_outdir in per_subject_eeg_outdirs.values())
    # Find per-subject unfold_out directories within EEG outdirs
    per_subject_unfold_out_result_dirs = {
        subject_id: os.path.join(subject_eeg_outdirs[0], "unfold_out")
        for subject_id, subject_eeg_outdirs in per_subject_eeg_outdirs.items()
    }

    # Iterate over subjects
    for subject_id in experiment.subject_ids:
        subject_unfold_out_result_dir = per_subject_unfold_out_result_dirs[subject_id]
        erp_noun_files = list_subject_files(
            subject_unfold_out_result_dir,
            subject_regex=subject_id,
            suffix=ERP_NOUN_FILE_SUFFIX
        )
        erp_fixation_files = list_subject_files(
            subject_unfold_out_result_dir,
            subject_regex=subject_id,
            suffix=ERP_FIXATION_FILE_SUFFIX
        )

        # Load noun data
        erp_noun_data = load_csv_list(
            erp_noun_files,
            progress_bar_description=f"Loading {len(erp_noun_files)} ERP noun files for subject <{subject_id}> ..."
        )
        erp_noun_data[CHANNEL_FIELD] = erp_noun_data[CHANNEL_FIELD].astype(str)

        # Load fixation data
        erp_fixation_data = load_csv_list(
            erp_fixation_files,
            progress_bar_description=f"Loading {len(erp_fixation_files)} ERP fixation files for subject <{subject_id}> ..."
        )
        erp_fixation_data[CHANNEL_FIELD] = erp_fixation_data[CHANNEL_FIELD].astype(str)

        # Merge channel coordinates into both dataframes
        erp_noun_data = erp_noun_data.merge(experiment.channel_coords, how="left", on=CHANNEL_FIELD)
        erp_fixation_data = erp_fixation_data.merge(experiment.channel_coords, how="left", on=CHANNEL_FIELD)

        # Annotate both dataframes for laterality and saggitality
        erp_noun_data, erp_fixation_data = map(annotate_laterality_and_saggitality, [erp_noun_data, erp_fixation_data])

    return experiment


if __name__ == "__main__":
    parser = argparse.ArgumentParser("Calculate and plot gaze statistics.")
    parser.add_argument('config', help='Path to config.yml file')
    args = parser.parse_args()
    main(os.path.abspath(args.config))
