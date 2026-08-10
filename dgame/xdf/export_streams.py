import argparse
import os

import pandas as pd
from scipy.io.wavfile import write as write_wav

from dgame.constants import BLOCK_IDS, ROUND_N
from dgame.xdf import (AUDIO_STREAM, EYETRACKER_STREAM,
                       SYNC_REFERENCE_BLOCK_COLUMN,
                       SYNC_REFERENCE_DRIFT_INTERCEPT_COLUMN,
                       SYNC_REFERENCE_DRIFT_SLOPE_COLUMN,
                       SYNC_REFERENCE_RAW_END_TIME_COLUMN,
                       SYNC_REFERENCE_RAW_START_TIME_COLUMN,
                       SYNC_REFERENCE_ROLE_COLUMN,
                       SYNC_REFERENCE_STREAM_COLUMN,
                       SYNC_REFERENCE_SYNCED_END_TIME_COLUMN,
                       SYNC_REFERENCE_SYNCED_START_TIME_COLUMN)
from dgame.xdf.utils import (extract_audio_stream_channels,
                             filter_streams_by_hostname)
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

        # Verify individual audio and time files
        subj_times_dir = os.path.join(experiment.times_outdir, subject_id)

        # One consolidated sync-reference file per subject (all blocks)
        sync_reference_file = os.path.join(subj_times_dir, f"{subject_id}_stream_sync.csv")
        assert_output_file_exists(sync_reference_file)

        for block in BLOCK_IDS:
            for role in experiment.participant_roles:
                role_audio_dir = experiment.get_role_outdir(experiment.audio_outdir, subject_id, role)
                audio_file = os.path.join(role_audio_dir, f"{subject_id}_{role}_{block}.wav")
                assert_output_file_exists(audio_file)


def get_eyetracker_streams_by_role(experiment, xdf: XDFFile) -> dict:
    """Get the eyetracker stream(s) present in an XDF recording, keyed by participant role.
    DGAME2 recordings contain a single (role-agnostic) eyetracker stream, keyed by None.
    DGAME3 recordings multiplex one eyetracker stream per dyad member into the same file
    (both members recorded on separate rigs), disambiguated by the hostname of each
    participant's assigned recording rig."""
    if experiment.dgame_version.n_participant_streams == 1:
        return {None: xdf.stream_by_name(EYETRACKER_STREAM)}

    streams_by_role = {}
    for role in experiment.participant_roles:
        hostname = experiment.get_rig_hostname(role)
        candidates = xdf.streams_by_name(EYETRACKER_STREAM)
        matches = filter_streams_by_hostname(candidates, hostname)
        if len(matches) != 1:
            raise ValueError(
                f"Expected exactly 1 <{EYETRACKER_STREAM}> stream for role <{role}> "
                f"(hostname <{hostname}>) in {xdf.path}, found {len(matches)}"
            )
        streams_by_role[role] = matches[0]
    return streams_by_role


def build_stream_sync_reference_rows(
        block: int,
        stream_label: str,
        stream: XDFStream,
        role: str = None,
    ) -> dict:
    """Assembles one sync-reference dict for a given block/stream, containing
    the stream's raw (own local clock) and LSL-synchronized start/end times, plus a
    linear clock-drift correction (synced = raw + drift_intercept + drift_slope * raw)
    fit from the stream's own clock-offset calibration data.
    `role` disambiguates DGAME3 recordings where more than one stream shares the same
    name (e.g. one eyetracker stream per dyad member); left as None for streams that
    aren't role-specific (e.g. the shared audio stream) or for DGAME2's single-role
    recordings, which have no such ambiguity to begin with.
    `stream` must come from an XDFFile loaded with synchronize_clocks=True."""
    drift_intercept, drift_slope = stream.fit_drift_correction()
    return {
        SYNC_REFERENCE_BLOCK_COLUMN: block,
        SYNC_REFERENCE_STREAM_COLUMN: stream_label,
        SYNC_REFERENCE_ROLE_COLUMN: role,
        SYNC_REFERENCE_RAW_START_TIME_COLUMN: round(stream.raw_start_time, ROUND_N),
        SYNC_REFERENCE_RAW_END_TIME_COLUMN: round(stream.raw_end_time, ROUND_N),
        SYNC_REFERENCE_SYNCED_START_TIME_COLUMN: round(stream.start_time, ROUND_N),
        SYNC_REFERENCE_SYNCED_END_TIME_COLUMN: round(stream.end_time, ROUND_N),
        SYNC_REFERENCE_DRIFT_INTERCEPT_COLUMN: round(drift_intercept, ROUND_N),
        # NB: drift_slope must not be rounded as it is a scalar, not a timestamp
        SYNC_REFERENCE_DRIFT_SLOPE_COLUMN: drift_slope,
    }


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
            xdf_file = experiment.get_xdf_file(subject_id, block, role=experiment.director_label)
            logger.info(f"Importing XDF file with clock synchronization: {xdf_file}")
            xdf_synced = XDFFile(
                xdf_file,
                synchronize_clocks=True,
                verbose=False,
            )

            # Extract audio stream channels to wav files, one per participant role
            # A single shared audio stream carries every role's channel,
            # for both DGAME2 and DGAME3 (only the per-role channel index differs)
            audio_stream_synced = xdf_synced.stream_by_name(AUDIO_STREAM)
            channel_samples, fs = extract_audio_stream_channels(audio_stream_synced)
            for role in experiment.participant_roles:
                samples = channel_samples[experiment.get_role_audio_channel(role)]
                role_audio_dir = experiment.get_role_outdir(experiment.audio_outdir, subject_id, role)
                os.makedirs(role_audio_dir, exist_ok=True)
                outwav = os.path.join(role_audio_dir, "_".join([subject_id, role, str(block)]) + ".wav")
                write_wav(outwav, int(fs), samples)

            # Record this block's raw and LSL-synchronized start/end times per stream,
            # so that later pipeline steps (which each only see one stream's
            # self-relative times) can recover the real inter-stream offset by looking
            # up a given stream's start time by name, instead of assuming all streams
            # started recording at the same instant.
            # Raw times come from the XDF footer chunk, which is unaffected by
            # clock synchronization, so a single synchronize_clocks=True load suffices;
            # no separate synchronize_clocks=False load of the same file is needed.
            # The audio stream is shared across roles, so it gets a single row without role specification;
            # the eyetracker stream gets one row per role for DGAME3 (disambiguated by hostname),
            # or a single row without role specification for DGAME2.
            subj_sync_reference_rows.append(
                build_stream_sync_reference_rows(block, AUDIO_STREAM, audio_stream_synced)
            )
            eyetracker_streams_by_role = get_eyetracker_streams_by_role(experiment, xdf_synced)
            for role, eyetracker_stream in eyetracker_streams_by_role.items():
                subj_sync_reference_rows.append(
                    build_stream_sync_reference_rows(block, EYETRACKER_STREAM, eyetracker_stream, role=role)
                )

        # Write one consolidated stream synchronization reference file per subject, covering all blocks
        sync_reference_csv = os.path.join(
            experiment.times_outdir,
            subject_id,
            f"{subject_id}_stream_sync.csv"
        )
        pd.DataFrame(subj_sync_reference_rows).to_csv(sync_reference_csv, index=False)

    # Validate outputs
    validate_outputs(experiment, experiment.subject_ids)

    return experiment


if __name__ == "__main__":
    parser = argparse.ArgumentParser("Export audio and eyetracking streams and create synchronization reference file.")
    parser.add_argument('config', help='Path to config.yml file')
    args = parser.parse_args()
    main(os.path.abspath(args.config))
