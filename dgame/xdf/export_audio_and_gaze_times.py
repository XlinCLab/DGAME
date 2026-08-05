import argparse
import os

import pandas as pd
from scipy.io.wavfile import write as write_wav

from dgame.constants import (BLOCK_IDS, DECKE_LABEL, DIRECTOR_LABEL,
                             PARTICIPANT_CONDITION_LABELS, ROUND_N)
from dgame.xdf import (AUDIO_STREAM, EYETRACKER_STREAM,
                       SYNC_REFERENCE_BLOCK_COLUMN,
                       SYNC_REFERENCE_DRIFT_INTERCEPT_COLUMN,
                       SYNC_REFERENCE_DRIFT_SLOPE_COLUMN,
                       SYNC_REFERENCE_RAW_END_TIME_COLUMN,
                       SYNC_REFERENCE_RAW_START_TIME_COLUMN,
                       SYNC_REFERENCE_STREAM_COLUMN,
                       SYNC_REFERENCE_SYNCED_END_TIME_COLUMN,
                       SYNC_REFERENCE_SYNCED_START_TIME_COLUMN)
from dgame.xdf.utils import extract_audio_stream_channels
from dgame.xdf.xdf_stream import XDFFile, XDFStream
from experiment.input_validation import (OutputValidationError,
                                         assert_output_file_exists)
from experiment.load_experiment import Experiment


def validate_outputs(experiment, subject_ids: list) -> None:
    """Validate audio output files."""
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

        # One consolidated sync-reference file per subject (all blocks)
        sync_reference_file = os.path.join(subj_times_dir, f"{subject_id}_sync_reference.csv")
        assert_output_file_exists(sync_reference_file)

        for block in BLOCK_IDS:
            for condition_label in PARTICIPANT_CONDITION_LABELS:
                audio_file = os.path.join(subject_audio_dir, f"{subject_id}_{condition_label}_{block}.wav")
                assert_output_file_exists(audio_file)


def build_stream_sync_reference_rows(
        block: int,
        audio_stream: XDFStream,
        eyetracker_stream: XDFStream,
    ) -> list[dict]:
    """Assembles one sync-reference row per stream for a given block, containing
    the stream's raw (own local clock) and LSL-synchronized start/end times, plus a
    linear clock-drift correction (synced = raw + drift_intercept + drift_slope * raw)
    fit from the stream's own clock-offset calibration data.
    Both streams must come from an XDFFile loaded with synchronize_clocks=True."""
    rows = []
    for stream_label, stream in (
        (AUDIO_STREAM, audio_stream),
        (EYETRACKER_STREAM, eyetracker_stream),
    ):
        drift_intercept, drift_slope = stream.fit_drift_correction()
        rows.append({
            SYNC_REFERENCE_BLOCK_COLUMN: block,
            SYNC_REFERENCE_STREAM_COLUMN: stream_label,
            SYNC_REFERENCE_RAW_START_TIME_COLUMN: round(stream.raw_start_time, ROUND_N),
            SYNC_REFERENCE_RAW_END_TIME_COLUMN: round(stream.raw_end_time, ROUND_N),
            SYNC_REFERENCE_SYNCED_START_TIME_COLUMN: round(stream.start_time, ROUND_N),
            SYNC_REFERENCE_SYNCED_END_TIME_COLUMN: round(stream.end_time, ROUND_N),
            SYNC_REFERENCE_DRIFT_INTERCEPT_COLUMN: round(drift_intercept, ROUND_N),
            # NB: drift_slope must not be rounded as it is a scalar, not a timestamp
            SYNC_REFERENCE_DRIFT_SLOPE_COLUMN: drift_slope,
        })
    return rows


def main(experiment: str | dict | Experiment) -> Experiment:
    # Initialize DGAME experiment from config
    if not isinstance(experiment, Experiment):
        from dgame.dgame import DGAME
        experiment = DGAME.from_input(experiment)
    logger = experiment.logger

    # Extract audio channels from XDF audio stream to wav files
    for subject_id in experiment.subject_ids:
        # Accumulated across all blocks and written once per subject, below
        subj_sync_reference_rows = []
        for block in BLOCK_IDS:
            xdf_file = experiment.get_xdf_file(subject_id, block, role=DIRECTOR_LABEL)
            logger.info(f"Importing XDF file with clock synchronization: {xdf_file}")
            xdf_synced = XDFFile(
                xdf_file,
                synchronize_clocks=True,
                verbose=False,
            )

            # Extract audio stream channels to wav files
            audio_stream_synced = xdf_synced.stream_by_name(AUDIO_STREAM)
            (director_samples, decke_samples), fs = extract_audio_stream_channels(audio_stream_synced)
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

            # Extract eyetracker stream
            eyetracker_stream_synced = xdf_synced.stream_by_name(EYETRACKER_STREAM)

            # Record this block's raw and LSL-synchronized start/end times per stream,
            # so that later pipeline steps (which each only see one stream's
            # self-relative times) can recover the real inter-stream offset by looking
            # up a given stream's start time by name, instead of assuming all streams
            # started recording at the same instant.
            # Raw times come from the XDF footer chunk, which is unaffected by
            # clock synchronization, so a single synchronize_clocks=True load suffices;
            # no separate synchronize_clocks=False load of the same file is needed.
            subj_sync_reference_rows.extend(build_stream_sync_reference_rows(
                block=block,
                audio_stream=audio_stream_synced,
                eyetracker_stream=eyetracker_stream_synced,
            ))

        # Write one consolidated sync-reference file per subject, covering all blocks
        sync_reference_csv = os.path.join(
            experiment.times_outdir,
            subject_id,
            f"{subject_id}_sync_reference.csv"
        )
        pd.DataFrame(subj_sync_reference_rows).to_csv(sync_reference_csv, index=False)

    # Validate outputs
    validate_outputs(experiment, experiment.subject_ids)

    return experiment


if __name__ == "__main__":
    parser = argparse.ArgumentParser("Export audio data and times of the gaze data.")
    parser.add_argument('config', help='Path to config.yml file')
    args = parser.parse_args()
    main(os.path.abspath(args.config))
