import argparse
import logging
import os

from pyxdf import load_xdf
from scipy.io.wavfile import write as write_wav

from dgame.constants import (AUDIO_STREAM, BLOCK_IDS, DECKE_LABEL,
                             DIRECTOR_LABEL, EYETRACKER_STREAM,
                             PARTICIPANT_CONDITION_LABELS, ROUND_N)
from experiment.input_validation import (OutputValidationError,
                                         assert_output_file_exists)
from experiment.load_experiment import Experiment
from utils.xdf_utils import (STREAM_TIMESTAMPS_LABEL,
                             extract_audio_stream_channels,
                             get_relative_times_from_stream, get_xdf_stream)

logger = logging.getLogger(__name__)


def validate_outputs(experiment, subject_ids: list) -> None:
    """Validate audio outputs from MATLAB script."""
    from dgame.dgame import validate_dgame_input
    experiment = validate_dgame_input(experiment)

    # Verify audio directories exist per subject and that expected files were created
    subject_audio_dirs_dict = experiment.get_subject_dirs_dict(experiment.audio_outdir)

    # Make sure audio and xdf directories are found for exactly the same subjects
    audio_subject_ids = sorted(list(subject_audio_dirs_dict.keys()))
    if audio_subject_ids != subject_ids:
        missing_audio = [subject_id for subject_id in subject_ids if subject_id not in subject_audio_dirs_dict]
        if len(missing_audio) > 0:
            raise OutputValidationError(f"Audio directory missing for following subjects: {', '.join(missing_audio)}")

    # Verify that expected audio and times files were created
    for subject_id, subject_audio_dirs in subject_audio_dirs_dict.items():
        # Verify that there is only one audio directory per subject
        try:
            assert len(subject_audio_dirs) == 1
        except AssertionError as exc:
            raise OutputValidationError(f">1 audio directory found for subject <{subject_id}>") from exc
        subject_audio_dir = subject_audio_dirs[0]

        # Verify individual audio and time files
        subj_times_dir = os.path.join(experiment.times_outdir, subject_id)
        for block in BLOCK_IDS:
            for condition_label in PARTICIPANT_CONDITION_LABELS:
                audio_file = os.path.join(subject_audio_dir, f"{subject_id}_{condition_label}_{block}.wav")
                assert_output_file_exists(audio_file)

                # time files per subject per block
                timestamp_file = os.path.join(subj_times_dir, f"{subject_id}_timestamps_max-min_{block}.txt")
                times_file = os.path.join(subj_times_dir, f"{subject_id}_times_{block}.txt")
                assert_output_file_exists(timestamp_file)
                assert_output_file_exists(times_file)


def main(experiment: str | dict | Experiment) -> Experiment:
    # Initialize DGAME experiment from config
    if not isinstance(experiment, Experiment):
        from dgame.dgame import DGAME
        experiment = DGAME.from_input(experiment)

    # Extract audio channels from XDF audio stream to wav files 
    for subject_id in experiment.subject_ids:
        subject_xdf_dir = os.path.join(experiment.xdf_indir, subject_id)
        for block in BLOCK_IDS:
            xdf_file = f"dgame{experiment.dgame_version}_{subject_id}_Director_{block}.xdf"
            xdf_file = os.path.join(subject_xdf_dir, "Director", xdf_file)
            xdf_data_with_clock_sync, _ = load_xdf(xdf_file, synchronize_clocks=True)

            # Extract audio stream channels to wav files
            audio_stream = get_xdf_stream(
                stream_label=AUDIO_STREAM,
                xdf_data=xdf_data_with_clock_sync,
            )
            (director_samples, decke_samples), fs = extract_audio_stream_channels(audio_stream)
            director_outwav = os.path.join(
                experiment.audio_outdir,
                subject_id,
                "_".join([subject_id, DIRECTOR_LABEL, str(block)]) + ".wav"
            )
            decke_outwav = os.path.join(
                experiment.audio_outdir,
                subject_id,
                "_".join([subject_id, DECKE_LABEL, str(block)]) + ".wav"
            )
            write_wav(director_outwav, int(fs), director_samples)
            write_wav(decke_outwav, int(fs), decke_samples)

            # Write eyetracker timestamps to CSV files
            # All timestamps as relative to first timestamp
            eyetracker_stream = get_xdf_stream(
                stream_label=EYETRACKER_STREAM,
                xdf_data=xdf_data_with_clock_sync,
            )
            relative_timestamps = get_relative_times_from_stream(
                eyetracker_stream,
                round_n=ROUND_N,  # TODO Ingmar's original MATLAB script rounded to 6 places instead of 7 (as default elsewhere), confirm whether this was intended
            )
            timestamp_csv = os.path.join(
                experiment.times_outdir,
                subject_id,
                "_".join([subject_id, "times", str(block)]) + ".txt"
            )
            with open(timestamp_csv, "w") as f:
                f.write("\n".join([str(t) for t in relative_timestamps]))

            # Get first and last timestamps rounded to ROUND_N decimal places
            # (NB: need to load XDF file without clock synchronization)
            xdf_data_no_clock_sync, _ = load_xdf(xdf_file, synchronize_clocks=False)
            eyetracker_stream = get_xdf_stream(
                stream_label=EYETRACKER_STREAM,
                xdf_data=xdf_data_no_clock_sync,
            )
            first_timestamp = round(eyetracker_stream[STREAM_TIMESTAMPS_LABEL][0], ROUND_N)
            last_timestamp = round(eyetracker_stream[STREAM_TIMESTAMPS_LABEL][-1], ROUND_N)
            max_min_timestamp_csv = os.path.join(
                experiment.times_outdir,
                subject_id,
                "_".join([subject_id, "timestamps", "max-min", str(block)]) + ".txt"
            )
            with open(max_min_timestamp_csv, "w") as f:
                f.write("\n".join([str(t) for t in [first_timestamp, last_timestamp]]))

    # Validate outputs
    validate_outputs(experiment, experiment.subject_ids)

    return experiment


if __name__ == "__main__":
    parser = argparse.ArgumentParser("Export audio data and times of the gaze data.")
    parser.add_argument('config', help='Path to config.yml file')
    args = parser.parse_args()
    main(os.path.abspath(args.config))
