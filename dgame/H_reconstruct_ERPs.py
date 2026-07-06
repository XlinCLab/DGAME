import argparse
import json
import os
from typing import Tuple

import numpy as np
import pandas as pd

from dgame.constants import (CONFLICT_LABEL, DET_POS_LABEL, FIXATION_LABEL,
                             NEXT_WORD_LABEL, NO_CONFLICT_LABEL,
                             NOUN_POS_LABEL, PREV_WORD_LABEL)
from experiment.input_validation import InputValidationError
from experiment.load_experiment import Experiment


def load_ufresult_struct(
        experiment: Experiment,
        ufresult_struct_file: str,
    ) -> Tuple[np.ndarray, np.ndarray, list[str], dict]:
    """
    Load `ufresult`-like output written by Unfold.jl (JLD2) and return (beta, times, chan_names).

    The file is created in Step G's Julia code as `*_ufresult_struct.jld2` and contains:
    - ufresult["beta"]: [channels x times x betas]
    - ufresult["times"]: time axis in seconds
    - ufresult["chanlocs"]: vector of dicts with at least {"labels": "<channel>"}
    - ufmeta: explicit coefficient metadata used to reconstruct MATLAB-like contrasts
    """
    jl = experiment.julia_interface
    jl.seval("using JLD2")
    jl.seval("using JSON3")
    path_lit = json.dumps(ufresult_struct_file)
    uf = jl.seval(
        f"""
        let
            d = JLD2.load({path_lit})
            uf = d["ufresult"]
            meta_json = JSON3.write(d["ufmeta"])
            beta = uf["beta"]
            times = uf["times"]
            chan_names = [String(c["labels"]) for c in uf["chanlocs"]]
            (beta, times, chan_names, meta_json)
        end
        """
    )
    beta = np.asarray(uf[0])
    times = np.asarray(uf[1]).reshape(-1)
    chan_names = [str(x) for x in list(uf[2])]
    uf_meta = json.loads(str(uf[3]))
    return beta, times, chan_names, uf_meta


def build_named_coef_vector(
        coef_order: list[str],
        assignments: dict[str, float],
        *,
        event_name: str,
        subject_id: str,
    ) -> np.ndarray:
    """Build a full coefficient vector by explicit name lookup."""
    v = np.zeros(len(coef_order), dtype=float)
    coef_index = {name: idx for idx, name in enumerate(coef_order)}
    qualified = {name: f"{event_name}::{name}" for name in assignments}
    missing = [name for name, qname in qualified.items() if qname not in coef_index]
    if missing:
        raise InputValidationError(
            f"Missing expected {event_name} coefficient(s) for subject {subject_id}: {', '.join(missing)}"
        )
    for name, value in assignments.items():
        v[coef_index[qualified[name]]] = float(value)
    return v


def load_events(events_file: str) -> pd.DataFrame:
    df = pd.read_csv(events_file)
    if "trial_time" in df.columns:
        df["trial_time"] = pd.to_numeric(df["trial_time"], errors="coerce")
    if "trial" in df.columns:
        df["trial"] = pd.to_numeric(df["trial"], errors="coerce")
    if "saccAmpl" in df.columns:
        df["saccAmpl"] = pd.to_numeric(df["saccAmpl"], errors="coerce")
    return df


def main(experiment: str | dict | Experiment) -> Experiment:
    # Initialize DGAME experiment from config
    if not isinstance(experiment, Experiment):
        from dgame.dgame import DGAME
        experiment = DGAME.from_input(experiment)
    logger = experiment.logger

    # Get list of subject IDs and their corresponding EEG directory paths
    subject_eeg_dirs_dict = experiment.get_subject_dirs_dict(experiment.eeg_outdir)
    subject_ids = list(subject_eeg_dirs_dict.keys())
    subject_eeg_dirs = list(subject_eeg_dirs_dict.values())
    assert all(len(subj_eeg_dir) == 1 for subj_eeg_dir in subject_eeg_dirs)
    subject_eeg_dirs = [subj_eeg_dir[0] for subj_eeg_dir in subject_eeg_dirs]

    # Python reimplementation of H_reconstruct_ERPs.m
    for subject_id, subject_eeg_dir in zip(subject_ids, subject_eeg_dirs):
        unfold_out_dir = os.path.join(subject_eeg_dir, "unfold_out")
        os.makedirs(unfold_out_dir, exist_ok=True)

        uf_struct_file = os.path.join(unfold_out_dir, f"{subject_id}_ufresult_struct.jld2")
        # Use the pre-unfold events CSV written by step G, not the original step F output.
        # Step G's update_fixation_events_df reclassifies conditionless fixations to `other_fixation`
        events_file = os.path.join(unfold_out_dir, f"{subject_id}_events_pre_unfold.csv")

        beta, times_s, chan_names, uf_meta = load_ufresult_struct(experiment, uf_struct_file)
        # Convert time axis to milliseconds for output CSVs
        times = times_s * 1000.0

        if beta.ndim == 2:
            # Julia exports coef(model) as [channels x (params * times)].
            # The coefficient layout is parameter-major within each channel, i.e. each
            # coefficient spans all time points before moving to the next coefficient.
            # Reshape to [channels x params x times] and transpose to [channels x times x params]
            # before reconstructing ERPs.
            n_channels, flat_len = beta.shape
            n_times = len(times_s)
            if flat_len % n_times != 0:
                raise InputValidationError(
                    f"Unexpected flattened beta length for subject {subject_id}: {flat_len} "
                    f"is not divisible by the number of time points ({n_times})"
                )
            n_params = flat_len // n_times
            beta = beta.reshape(n_channels, n_params, n_times).transpose(0, 2, 1)
        elif beta.ndim == 3:
            n_channels, n_times, n_params = beta.shape
        else:
            raise InputValidationError(f"Unexpected beta array shape for subject {subject_id}: {beta.shape}")

        if n_channels != len(chan_names):
            raise InputValidationError(
                f"Channel count mismatch for subject {subject_id}: beta has {n_channels} channels "
                f"but ufresult has {len(chan_names)} chan labels"
            )

        if "coef_order" not in uf_meta or "coefnames_by_event" not in uf_meta or "event_order" not in uf_meta:
            raise InputValidationError(
                f"Missing coefficient metadata in Julia output for subject {subject_id}"
            )
        coef_order = [str(name) for name in uf_meta["coef_order"]]
        coefnames_by_event = {
            str(event): [str(name) for name in names]
            for event, names in uf_meta["coefnames_by_event"].items()
        }
        event_order = [str(name) for name in uf_meta["event_order"]]
        if len(coef_order) != n_params:
            raise InputValidationError(
                f"Coefficient metadata mismatch for subject {subject_id}: "
                f"metadata has {len(coef_order)} coefficients but beta has {n_params}"
            )
        expected_events = {PREV_WORD_LABEL, NEXT_WORD_LABEL, FIXATION_LABEL, DET_POS_LABEL, NOUN_POS_LABEL}
        missing_events = expected_events.difference(coefnames_by_event)
        if missing_events:
            raise InputValidationError(
                f"Missing expected event block(s) in Julia metadata for subject {subject_id}: "
                f"{', '.join(sorted(missing_events))}"
            )
        if event_order != [PREV_WORD_LABEL, NEXT_WORD_LABEL, FIXATION_LABEL, DET_POS_LABEL, NOUN_POS_LABEL]:
            logger.info(
                f"Subject {subject_id}: Julia event order is {event_order}; "
                "Step H will reconstruct by explicit coefficient name."
            )

        events = load_events(events_file)
        mean_sacc = np.nanmean(events.loc[events["type"] == FIXATION_LABEL, "saccAmpl"].to_numpy())
        if np.isnan(mean_sacc):
            mean_sacc = 0.0

        con_idx = (events["type"] == NOUN_POS_LABEL) & (events["condition"] == CONFLICT_LABEL)
        ncon_idx = (events["type"] == NOUN_POS_LABEL) & (events["condition"] == NO_CONFLICT_LABEL)
        tg_con = np.unique(events.loc[con_idx, "trial_time"].dropna().to_numpy())
        tg_nocon = np.unique(events.loc[ncon_idx, "trial_time"].dropna().to_numpy())
        trial_con = np.nanmean(events.loc[con_idx, "trial"].to_numpy())
        trial_nocon = np.nanmean(events.loc[ncon_idx, "trial"].to_numpy())

        # Noun ERPs
        all_times_n = np.concatenate([tg_con, tg_nocon]) if tg_con.size + tg_nocon.size > 0 else np.array([])
        is_con_n = np.concatenate([np.ones(len(tg_con), dtype=bool), np.zeros(len(tg_nocon), dtype=bool)])

        for ch_idx, ch_name in enumerate(chan_names):
            beta_mat = beta[ch_idx, :, :]  # [times x params]
            if all_times_n.size == 0:
                continue

            erp_n = []
            erp_records = []
            for tval, is_con in zip(all_times_n, is_con_n):
                conflict = 1.0 if is_con else 0.0
                trial_mean = float(trial_con if is_con else trial_nocon)

                # MATLAB reconstruction logic for noun ERPs:
                #   - intercept is always on
                #   - the condition indicator is on only for conflict nouns
                #   - trial_time uses the actual noun trial time on that row
                #   - trial uses the mean trial number for that condition
                #   - the interaction equals trial_time for conflict nouns, else zero
                # This mirrors the MATLAB Step H vector construction, but assigns each
                # coefficient explicitly by event-qualified Julia coefficient name.
                assignments = {
                    "(Intercept)": 1.0,
                    "condition: conflict": conflict,
                    "trial_time": float(tval),
                    "trial": trial_mean,
                    "condition: conflict & trial_time": float(tval) if conflict else 0.0,
                }
                v = build_named_coef_vector(
                    coef_order,
                    assignments,
                    event_name=NOUN_POS_LABEL,
                    subject_id=subject_id,
                )

                erp = beta_mat @ v  # [times]
                erp_n.append(erp)
                erp_records.append((CONFLICT_LABEL if is_con else NO_CONFLICT_LABEL, tval))

            erp_n = np.column_stack(erp_n)
            n_rows = erp_n.shape[0] * erp_n.shape[1]
            T_time = np.tile(times, erp_n.shape[1])
            T_event = np.repeat("noun", n_rows)
            T_cond = np.repeat([r[0] for r in erp_records], erp_n.shape[0])
            T_mtf = np.repeat([r[1] for r in erp_records], erp_n.shape[0])
            T_chan = np.repeat(ch_name, n_rows)
            # Betas are already in µV: step G converts MNE Volts → µV before passing
            # data to Unfold.jl, so the fitted coefficients carry µV units throughout.
            T_data = erp_n.reshape(-1, order="F")
            T_subj = np.repeat(subject_id, n_rows)

            noun_table = pd.DataFrame(
                {
                    "time": T_time,
                    "event": T_event,
                    "condition": T_cond,
                    "mean_target_fixation": T_mtf,
                    "channel": T_chan,
                    "data": T_data,
                    "subject": T_subj,
                }
            )
            noun_table.to_csv(
                os.path.join(unfold_out_dir, f"{subject_id}_{ch_name}_unfold_N.csv"),
                index=False,
            )

        # Fixation ERPs
        fix_times = np.unique(np.round(events.loc[events["type"] == FIXATION_LABEL, "trial_time"].dropna(), 1))
        fix_labels = ["elsewhere", "other", "target"]
        con_vals = [True, False]

        for ch_idx, ch_name in enumerate(chan_names):
            beta_mat = beta[ch_idx, :, :]
            erp_list = []
            erp_records = []
            for ft in fix_times:
                for cv in con_vals:
                    for fl in fix_labels:
                        conflict = 1.0 if cv else 0.0
                        is_other = 1.0 if fl == "other" else 0.0
                        is_target = 1.0 if fl == "target" else 0.0

                        # MATLAB reconstruction logic for fixation ERPs:
                        #   - intercept always on
                        #   - condition, fix_at, and trial_time terms follow the observed
                        #     fixation condition/label/time
                        #   - interaction terms are the corresponding products
                        #   - spline terms are set to the mean saccadic amplitude
                        # We keep the logic explicit and name-based so it remains robust
                        # to any future Unfold coefficient reordering.
                        assignments = {
                            "(Intercept)": 1.0,
                            "condition: conflict": conflict,
                            "fix_at: other": is_other,
                            "fix_at: target": is_target,
                            "trial_time": ft,
                            "condition: conflict & fix_at: other": conflict * is_other,
                            "condition: conflict & fix_at: target": conflict * is_target,
                            "condition: conflict & trial_time": conflict * ft,
                            "fix_at: other & trial_time": is_other * ft,
                            "fix_at: target & trial_time": is_target * ft,
                            "condition: conflict & fix_at: other & trial_time": conflict * is_other * ft,
                            "condition: conflict & fix_at: target & trial_time": conflict * is_target * ft,
                        }
                        for coefname in coefnames_by_event[FIXATION_LABEL]:
                            if coefname.startswith("spl("):
                                assignments[coefname] = mean_sacc
                        v = build_named_coef_vector(
                            coef_order,
                            assignments,
                            event_name=FIXATION_LABEL,
                            subject_id=subject_id,
                        )

                        erp = beta_mat @ v  # [times]
                        erp_list.append(erp)
                        erp_records.append((CONFLICT_LABEL if cv else NO_CONFLICT_LABEL, fl, float(ft)))

            if not erp_list:
                continue
            erp_f = np.column_stack(erp_list)
            n_rows = erp_f.shape[0] * erp_f.shape[1]
            T_time = np.tile(times, erp_f.shape[1])
            T_event = np.repeat(FIXATION_LABEL, n_rows)
            T_cond = np.repeat([r[0] for r in erp_records], erp_f.shape[0])
            T_fixat = np.repeat([r[1] for r in erp_records], erp_f.shape[0])
            T_trtime = np.repeat([r[2] for r in erp_records], erp_f.shape[0])
            T_chan = np.repeat(ch_name, n_rows)
            # Betas are already in µV: step G converts MNE Volts → µV before passing
            # data to Unfold.jl, so the fitted coefficients carry µV units throughout.
            T_data = erp_f.reshape(-1, order="F")
            T_subj = np.repeat(subject_id, n_rows)

            fix_table = pd.DataFrame(
                {
                    "time": T_time,
                    "event": T_event,
                    "condition": T_cond,
                    "fix_at": T_fixat,
                    "trial_time": T_trtime,
                    "channel": T_chan,
                    "data": T_data,
                    "subject": T_subj,
                }
            )
            fix_table.to_csv(
                os.path.join(unfold_out_dir, f"{subject_id}_{ch_name}_unfold_FIX.csv"),
                index=False,
            )

    return experiment


if __name__ == "__main__":
    parser = argparse.ArgumentParser("Reconstruct ERPs using the GLM beta estimates and export to CSV.")
    parser.add_argument('config', help='Path to config.yml file')
    args = parser.parse_args()
    main(os.path.abspath(args.config))
