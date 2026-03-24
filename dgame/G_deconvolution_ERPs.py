import argparse
import copy
import os
from typing import Any

import numpy as np
from juliacall import Main as jl
from juliacall import Pkg as jlPkg
from scipy.io import loadmat, savemat

from dgame.constants import (CONFLICT_LABEL, NO_CONFLICT_LABEL, STEP_F_KEY,
                             TRIAL_TIME_OFFSET)
from experiment.input_validation import InputValidationError
from experiment.load_experiment import Experiment


def _to_str(value: Any) -> str | None:
    if value is None:
        return None
    if isinstance(value, bytes):
        return value.decode("utf-8")
    if isinstance(value, (list, tuple)):
        if len(value) == 0:
            return None
        return _to_str(value[0])
    if isinstance(value, np.ndarray):
        if value.size == 0:
            return None
        if value.size == 1:
            return _to_str(value.item())
        if value.dtype.kind in {"U", "S"}:
            return "".join([_to_str(v) or "" for v in value.tolist()])
        return _to_str(value.flat[0])
    if isinstance(value, str):
        return value
    return str(value)


def _to_float(value: Any) -> float | None:
    if value is None:
        return None
    if isinstance(value, np.ndarray):
        if value.size == 0:
            return None
        value = value.item() if value.size == 1 else value.flat[0]
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def _get_field(event: Any, name: str, default: Any = None) -> Any:
    if isinstance(event, dict):
        return event.get(name, default)
    return getattr(event, name, default)


def _set_field(event: Any, name: str, value: Any) -> None:
    if isinstance(event, dict):
        event[name] = value
    else:
        setattr(event, name, value)


def _is_type(event: Any, type_value: str) -> bool:
    return (_to_str(_get_field(event, "type")) or "") == type_value


def _as_event_list(events: Any) -> list[Any]:
    if events is None:
        return []
    if isinstance(events, np.ndarray):
        return list(events.ravel())
    return [events]


def _update_fixation_events(events: list[Any]) -> None:
    # setup the fixations for analysis
    # loop through events and find fixations to add condition information retrieved
    # from the entry of the noun of the corresponding trial
    for ev_idx, ev in enumerate(events):
        if not _is_type(ev, "fixation"):
            continue

        # summarize non-target fixations
        fix_at = _to_str(_get_field(ev, "fix_at")) or ""
        if fix_at not in {"target", "other"}:
            _set_field(ev, "fix_at", "elsewhere")

        # if it has trial time value smaller 0, i.e. occurs before noun onset
        trial_time = _to_float(_get_field(ev, "trial_time"))
        if trial_time is not None and trial_time < 0 and trial_time >= -TRIAL_TIME_OFFSET:
            # go forward a few events and see if you find a noun, stop at the first one
            # (200 is arbitrary so we are sure to get all nouns in all trials)
            for find_idx in range(1, 201):
                cand_idx = ev_idx + find_idx
                if cand_idx >= len(events):
                    break
                cand = events[cand_idx]
                if _is_type(cand, "N"):
                    _set_field(ev, "condition", _get_field(cand, "condition"))
                    _set_field(ev, "set", _get_field(cand, "set"))
                    _set_field(ev, "pattern", _get_field(cand, "pattern"))
                    break
        elif trial_time is not None and trial_time > 0 and trial_time <= TRIAL_TIME_OFFSET:
            # do the same for fixations after noun onset
            # go back a few events and see if you find a noun, stop at the first one
            # (200 is arbitrary so we are sure to get all nouns in all trials)
            for find_idx in range(1, 201):
                cand_idx = ev_idx - find_idx
                if cand_idx < 0:
                    break
                cand = events[cand_idx]
                if _is_type(cand, "N"):
                    _set_field(ev, "condition", _get_field(cand, "condition"))
                    _set_field(ev, "set", _get_field(cand, "set"))
                    _set_field(ev, "pattern", _get_field(cand, "pattern"))
                    break

        condition = _to_str(_get_field(ev, "condition")) or ""
        if condition not in {CONFLICT_LABEL, NO_CONFLICT_LABEL}:
            # set all fixation without a condition to a different type
            _set_field(ev, "type", "other_fixation")
            _set_field(ev, "condition_old", condition or "NA")
            _set_field(ev, "condition", "NA")

    # ensure fixation fix_at labeling is consistent
    for ev in events:
        if not _is_type(ev, "fixation"):
            continue
        fix_at = _to_str(_get_field(ev, "fix_at")) or ""
        if fix_at not in {"target", "other"}:
            _set_field(ev, "fix_at", "elsewhere")


def _load_eeglab_set(set_path: str) -> Any:
    try:
        mat = loadmat(set_path, struct_as_record=False, squeeze_me=True)
    except NotImplementedError as exc:
        raise RuntimeError(
            f"Unable to load EEGLAB .set file {set_path}. "
            "If this file was saved as MATLAB v7.3 (HDF5), "
            "resave as v7.0/v7.2 or install h5py and load via HDF5." 
        ) from exc

    eeg = mat.get("EEG")
    if eeg is not None:
        return eeg

    # Some EEGLAB .set files store fields at top-level (no EEG struct)
    # Treat the full dict (minus MATLAB metadata) as the EEG struct
    eeg = {
        key: value
        for key, value in mat.items()
        if key not in {"__header__", "__version__", "__globals__"}
    }
    if not eeg:
        raise RuntimeError(f"No EEG struct found in {set_path}")
    return eeg


def _save_eeglab_set(set_path: str, eeg: Any) -> None:
    savemat(set_path, {"EEG": eeg}, do_compression=True)


def main(experiment: str | dict | Experiment) -> Experiment:
    # Initialize DGAME experiment from config
    if not isinstance(experiment, Experiment):
        from dgame.dgame import DGAME
        experiment = DGAME.from_input(experiment)
    logger = experiment.logger

    # Get experiment subject IDs and their corresponding EEG directory paths
    subject_eeg_dirs_dict = experiment.get_subject_dirs_dict(experiment.eeg_outdir)

    # Run pre-unfold processing
    for subject_id, subject_dirs in subject_eeg_dirs_dict.items():
        subject_eeg_dir = subject_dirs[0]
        cleaned_set = f"{subject_id}_director_cleaned.set"
        cleaned_set_path = os.path.join(subject_eeg_dir, cleaned_set)
        if not os.path.exists(cleaned_set_path):
            raise InputValidationError(f"Missing cleaned EEG set (from step {STEP_F_KEY}): {cleaned_set_path}")
        logger.info(f"Loading cleaned EEG set for subject {subject_id}: {cleaned_set_path}")
        eeg = _load_eeglab_set(cleaned_set_path)

        # Backup event data as attribute of new EEG object 
        event_old = copy.deepcopy(_get_field(eeg, "event"))
        _set_field(eeg, "event_old", event_old)

        events = _as_event_list(_get_field(eeg, "event"))
        if not events:
            logger.warning(f"No events found in EEG set for subject {subject_id}")
            continue

        _update_fixation_events(events)

        pre_unfold_set = os.path.join(subject_eeg_dir, f"{subject_id}_director_cleaned_pre_unfold.set")
        logger.info(f"Saving pre-unfold EEG set for subject {subject_id}: {pre_unfold_set}")
        _save_eeglab_set(pre_unfold_set, eeg)

        # output directory for unfold results
        outpath = os.path.join(subject_eeg_dir, "unfold_out")
        os.makedirs(outpath, exist_ok=True)

    # Run Unfold analysis in Julia
    logger.info("Running unfold analysis in Julia...")
    unfold_julia_script = os.path.join(experiment.julia_script_dir, "unfold_step_g.jl")
    jlPkg.activate(experiment.julia_script_dir)
    jlPkg.instantiate()
    jl.seval(f'include("{unfold_julia_script}")')
    for subject_id, subject_dirs in subject_eeg_dirs_dict.items():
        subject_eeg_dir = subject_dirs[0]
        pre_unfold_set = os.path.join(subject_eeg_dir, f"{subject_id}_director_cleaned_pre_unfold.set")
        outpath = os.path.join(subject_eeg_dir, "unfold_out")
        logger.info(f"Running Unfold.jl for subject {subject_id}...")
        jl.run_unfold_step_g(pre_unfold_set, outpath, subject_id)

    return experiment


if __name__ == "__main__":
    parser = argparse.ArgumentParser("Compute rERPs.")
    parser.add_argument('config', help='Path to config.yml file')
    args = parser.parse_args()
    main(os.path.abspath(args.config))
