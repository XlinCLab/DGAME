import os

from dgame.constants import BLOCK_IDS


def retrieve_and_validate_inputs(experiment) -> tuple[list, list]:
    """Validate that all required input files exist and return flattened lists of subject IDs and per-subject input directories."""
    from dgame.dgame import validate_dgame_input
    experiment = validate_dgame_input(experiment)

    # Get XDF directory paths per subject
    subject_xdf_dirs_dict = experiment.get_subject_dirs_dict(experiment.xdf_indir)
    subject_xdf_dir_list = []
    subject_ids = []
    for subject_id, subject_xdf_dirs in subject_xdf_dirs_dict.items():
        # Verify that there is only one xdf directory per subject
        try:
            assert len(subject_xdf_dirs) == 1
        except AssertionError as exc:
            raise ValueError(f">1 xdf directory found for subject <{subject_id}>") from exc
        subject_xdf_dir = subject_xdf_dirs[0]
        subject_xdf_dir_list.append(subject_xdf_dir)
        subject_ids.append(subject_id)
        
        # Verify that the xdf directory contains all required files
        subject_xdf_director_dir = os.path.join(subject_xdf_dir, "Director")
        for block in BLOCK_IDS:
            xdf_file = os.path.join(subject_xdf_director_dir, f"dgame2_{subject_id}_Director_{str(block)}.xdf")
            try:
                assert os.path.exists(xdf_file)
            except AssertionError as exc:
                raise FileNotFoundError(f"Input file {xdf_file} not found")
    subject_ids.sort()

    # Create audio directories per subject
    for subject_id in subject_ids:
        subject_audio_dir = os.path.join(experiment.audio_indir, subject_id)
        os.makedirs(subject_audio_dir, exist_ok=True)
    
    return subject_ids, subject_xdf_dir_list
