import argparse
import logging
import os

from pyxdf import load_xdf
from scipy.io.wavfile import write as write_wav

from dgame.constants import AUDIO_STREAM, BLOCK_IDS, EYETRACKER_STREAM, ROUND_N
from experiment.input_validation import (OutputValidationError,
                                         assert_output_file_exists)
from experiment.load_experiment import Experiment
from utils.xdf_utils import (STREAM_TIMESTAMPS_LABEL,
                             extract_audio_stream_channels,
                             fill_stream_gaps, filter_streams_by_hostname,
                             get_relative_times_from_stream, get_xdf_stream,
                             get_xdf_streams_by_name)

# DGAME2 recordings only ever contain a single eyetracker stream (only one dyad
# member was ET/EEG-equipped); DGAME3 recordings contain one eyetracker stream
# per dyad member, disambiguated by the hostname of their fixedly assigned rig
SINGLE_PARTICIPANT_DGAME_VERSION = "2"


def get_role_outdir(base_dir: str, subject_id: str, role: str, dgame_version: str) -> str:
    """Get the (optionally role-specific) output directory for a given subject and role."""
    if dgame_version == SINGLE_PARTICIPANT_DGAME_VERSION:
        return os.path.join(base_dir, subject_id)
    return os.path.join(base_dir, subject_id, role)


def get_eyetracker_streams_by_role(experiment, xdf_data: list, xdf_file: str) -> dict:
    """Get the eyetracker stream(s) present in an XDF recording, keyed by participant role.
    For dgame2, there is a single (role-agnostic) eyetracker stream. For dgame3, there is one
    eyetracker stream per dyad member, disambiguated via the hostname of their assigned rig."""
    if experiment.dgame_version == SINGLE_PARTICIPANT_DGAME_VERSION:
        return {None: get_xdf_stream(stream_label=EYETRACKER_STREAM, xdf_data=xdf_data)}

    streams_by_role = {}
    for role in experiment.participant_roles:
        hostname = experiment.get_rig_hostname(role)
        candidates = get_xdf_streams_by_name(EYETRACKER_STREAM, xdf_data=xdf_data)
        matches = filter_streams_by_hostname(candidates, hostname)
        if len(matches) != 1:
            raise ValueError(
                f"Expected exactly 1 <{EYETRACKER_STREAM}> stream for role <{role}> "
                f"(hostname <{hostname}>) in {xdf_file}, found {len(matches)}"
            )
        streams_by_role[role] = matches[0]
    return streams_by_role


def get_block_output_files(experiment, subject_id: str, block: int) -> list[str]:
    """Get the full list of expected output files (audio + times) for a given subject/block."""
    output_files = []
    for role in experiment.participant_roles:
        role_audio_dir = get_role_outdir(experiment.audio_outdir, subject_id, role, experiment.dgame_version)
        output_files.append(os.path.join(role_audio_dir, f"{subject_id}_{role}_{block}.wav"))

        role_times_dir = get_role_outdir(experiment.times_outdir, subject_id, role, experiment.dgame_version)
        output_files.append(os.path.join(role_times_dir, f"{subject_id}_times_{block}.txt"))
        output_files.append(os.path.join(role_times_dir, f"{subject_id}_timestamps_max-min_{block}.txt"))
    return output_files


def prompt_skip_existing_blocks(already_done: list[tuple]) -> bool:
    """Ask the user once whether to skip all blocks that already have complete output, or
    reprocess them anyway. Returns True to skip, False to reprocess."""
    block_labels = ", ".join(f"{subject_id}/{block}" for subject_id, block, _ in already_done)
    print(f"{len(already_done)} block(s) already have complete output, skipping: {block_labels}")
    response = input("Skip (Y) or run again (R)? [Y/r]: ").strip().lower()
    return response != "r"


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
        for block in BLOCK_IDS:
            for role in experiment.participant_roles:
                role_audio_dir = get_role_outdir(experiment.audio_outdir, subject_id, role, experiment.dgame_version)
                audio_file = os.path.join(role_audio_dir, f"{subject_id}_{role}_{block}.wav")
                assert_output_file_exists(audio_file)

                # time files per subject per block (per role, for dgame3; shared for dgame2)
                role_times_dir = get_role_outdir(experiment.times_outdir, subject_id, role, experiment.dgame_version)
                timestamp_file = os.path.join(role_times_dir, f"{subject_id}_timestamps_max-min_{block}.txt")
                times_file = os.path.join(role_times_dir, f"{subject_id}_times_{block}.txt")
                assert_output_file_exists(timestamp_file)
                assert_output_file_exists(times_file)


def main(experiment: str | dict | Experiment) -> Experiment:
    # Initialize DGAME experiment from config
    if not isinstance(experiment, Experiment):
        from dgame.dgame import DGAME
        experiment = DGAME.from_input(experiment)
    logger = experiment.logger

    # Find subject/blocks whose output already fully exists, and ask once whether to skip
    # them (avoids redoing expensive XDF re-loads for already-completed blocks on a rerun)
    already_done = [
        (subject_id, block, os.path.join(experiment.xdf_indir, subject_id, "Director",
                                         f"dgame{experiment.dgame_version}_{subject_id}_Director_{block}.xdf"))
        for subject_id in experiment.subject_ids
        for block in BLOCK_IDS
        if all(os.path.exists(f) for f in get_block_output_files(experiment, subject_id, block))
    ]
    skip_existing = prompt_skip_existing_blocks(already_done) if already_done else False
    already_done_blocks = {(subject_id, block) for subject_id, block, _ in already_done}

    # Extract audio channels from XDF audio stream to wav files
    failed_xdf_files = []
    for subject_id in experiment.subject_ids:
        subject_xdf_dir = os.path.join(experiment.xdf_indir, subject_id)
        for block in BLOCK_IDS:
            xdf_file = f"dgame{experiment.dgame_version}_{subject_id}_Director_{block}.xdf"
            xdf_file = os.path.join(subject_xdf_dir, "Director", xdf_file)

            if skip_existing and (subject_id, block) in already_done_blocks:
                logger.info(f"Output already exists for {xdf_file}, skipping")
                continue

            try:
                logger.info(f"Importing XDF file with clock synchronization: {xdf_file}")
                xdf_data_with_clock_sync, _ = load_xdf(
                    xdf_file,
                    synchronize_clocks=True,
                    verbose=False,
                )

                # Extract audio stream channels to wav files
                # NB: gaps between jitter-removal segments (see pyxdf's "Segments and
                # clock-segments differ" warning) are filled with zero-valued (silent)
                # samples first, so the exported audio stays aligned with other streams
                # instead of silently drifting shorter than the real recording duration
                audio_stream = get_xdf_stream(
                    stream_label=AUDIO_STREAM,
                    xdf_data=xdf_data_with_clock_sync,
                )
                audio_stream = fill_stream_gaps(audio_stream)
                channel_samples, fs = extract_audio_stream_channels(audio_stream)
                if len(channel_samples) == 0:
                    raise ValueError(f"Audio stream in {xdf_file} contains 0 samples (recording is empty)")
                for role in experiment.participant_roles:
                    samples = channel_samples[experiment.get_role_audio_channel(role)]
                    role_audio_dir = get_role_outdir(experiment.audio_outdir, subject_id, role, experiment.dgame_version)
                    os.makedirs(role_audio_dir, exist_ok=True)
                    outwav = os.path.join(role_audio_dir, "_".join([subject_id, role, str(block)]) + ".wav")
                    write_wav(outwav, int(fs), samples)

                # Write eyetracker timestamps to CSV files
                # All timestamps as relative to first timestamp
                eyetracker_streams_by_role = get_eyetracker_streams_by_role(experiment, xdf_data_with_clock_sync, xdf_file)
                for role, eyetracker_stream in eyetracker_streams_by_role.items():
                    relative_timestamps = get_relative_times_from_stream(
                        eyetracker_stream,
                        round_n=ROUND_N,
                    )
                    role_times_dir = get_role_outdir(experiment.times_outdir, subject_id, role, experiment.dgame_version)
                    os.makedirs(role_times_dir, exist_ok=True)
                    timestamp_csv = os.path.join(
                        role_times_dir,
                        "_".join([subject_id, "times", str(block)]) + ".txt"
                    )
                    with open(timestamp_csv, "w") as f:
                        f.write("\n".join([str(t) for t in relative_timestamps]))

                # Get first and last timestamps rounded to ROUND_N decimal places
                # (NB: need to load XDF file without clock synchronization)
                logger.info(f"Importing XDF file without clock synchronization: {xdf_file}")
                xdf_data_no_clock_sync, _ = load_xdf(
                    xdf_file,
                    synchronize_clocks=False,
                    verbose=False,
                )
                eyetracker_streams_by_role = get_eyetracker_streams_by_role(experiment, xdf_data_no_clock_sync, xdf_file)
                for role, eyetracker_stream in eyetracker_streams_by_role.items():
                    first_timestamp = round(eyetracker_stream[STREAM_TIMESTAMPS_LABEL][0], ROUND_N)
                    last_timestamp = round(eyetracker_stream[STREAM_TIMESTAMPS_LABEL][-1], ROUND_N)
                    role_times_dir = get_role_outdir(experiment.times_outdir, subject_id, role, experiment.dgame_version)
                    max_min_timestamp_csv = os.path.join(
                        role_times_dir,
                        "_".join([subject_id, "timestamps", "max-min", str(block)]) + ".txt"
                    )
                    with open(max_min_timestamp_csv, "w") as f:
                        f.write("\n".join([str(t) for t in [first_timestamp, last_timestamp]]))
            except Exception as exc:
                # A single unrecordable/corrupted block (e.g. an empty audio stream) must not
                # abort export for every other subject/block; log it and keep going -
                # validate_outputs() below will surface exactly which outputs are still missing
                logger.error(f"Failed to export audio/times for {xdf_file}: {exc}")
                failed_xdf_files.append(xdf_file)
                continue

    if failed_xdf_files:
        logger.warning(
            f"Failed to process {len(failed_xdf_files)} XDF file(s), skipped: "
            + ", ".join(failed_xdf_files)
        )

    # Validate outputs
    validate_outputs(experiment, experiment.subject_ids)

    return experiment


if __name__ == "__main__":
    parser = argparse.ArgumentParser("Export audio data and times of the gaze data.")
    parser.add_argument('config', help='Path to config.yml file')
    args = parser.parse_args()
    main(os.path.abspath(args.config))
