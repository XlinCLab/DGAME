import argparse
import os
from typing import Tuple

import numpy as np
import pandas as pd

from dgame.constants import STEP_H_KEY
from experiment.load_experiment import Experiment


def _load_ufresult_struct_from_julia(
    experiment: Experiment,
    ufresult_struct_file: str,
) -> Tuple[np.ndarray, np.ndarray, list[str]]:
    """
    Load `ufresult`-like output written by Unfold.jl (JLD2) and return (beta, times, chan_names).

    The file is created in Step G's Julia code as `*_ufresult_struct.jld2` and contains:
    - ufresult["beta"]: [channels x times x betas]
    - ufresult["times"]: time axis in seconds
    - ufresult["chanlocs"]: vector of dicts with at least {"labels": "<channel>"}
    """
    jl = experiment.julia_interface
    jl.seval("using JLD2")
    uf = jl.seval(
        f"""
        let
            d = JLD2.load({ufresult_struct_file!r})
            uf = d["ufresult"]
            beta = uf["beta"]
            times = uf["times"]
            chan_names = [String(c["labels"]) for c in uf["chanlocs"]]
            (beta, times, chan_names)
        end
        """
    )
    beta = np.asarray(uf[0])
    times = np.asarray(uf[1]).reshape(-1)
    chan_names = [str(x) for x in list(uf[2])]
    return beta, times, chan_names


def _load_events(events_file: str) -> pd.DataFrame:
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
        events_file = os.path.join(subject_eeg_dir, f"{subject_id}_director_events.csv")

        beta, times_s, chan_names = _load_ufresult_struct_from_julia(experiment, uf_struct_file)
        # Convert time axis to milliseconds for output CSVs
        times = times_s * 1000.0

        if beta.ndim != 3:
            raise RuntimeError(f"Unexpected beta array shape for subject {subject_id}: {beta.shape}")

        n_channels, n_times, n_params = beta.shape
        if n_channels != len(chan_names):
            raise RuntimeError(
                f"Channel count mismatch for subject {subject_id}: beta has {n_channels} channels "
                f"but ufresult has {len(chan_names)} chan labels"
            )

        # Unfold design in Step G Julia script:
        # prev: 1, next: 1, fixation: 12 + k_spline, D: 3, N: 5  => total = 22 + k_spline
        k_spline = n_params - 22
        if k_spline < 0:
            raise RuntimeError(
                f"Unexpected number of parameters for subject {subject_id}: {n_params} (expected >= 22)"
            )

        events = _load_events(events_file)
        mean_sacc = np.nanmean(events.loc[events["type"] == "fixation", "saccAmpl"].to_numpy())
        if np.isnan(mean_sacc):
            mean_sacc = 0.0

        con_idx = (events["type"] == "N") & (events["condition"] == "conflict")
        ncon_idx = (events["type"] == "N") & (events["condition"] == "no_conflict")
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
            meta = []
            for tval, is_con in zip(all_times_n, is_con_n):
                # Construct a covariate vector aligned to the Unfold coefficient order.
                # We only populate the N (noun) block (last 5 params) and leave others at 0.
                v = np.zeros(n_params, dtype=float)
                base = n_params - 5
                conflict = 1.0 if is_con else 0.0
                trial_mean = float(trial_con if is_con else trial_nocon)

                # N formula in Step G Julia script: 0 ~ 1 + condition * trial_time + trial
                v[base + 0] = 1.0
                v[base + 1] = conflict
                v[base + 2] = float(tval)
                v[base + 3] = conflict * float(tval)
                v[base + 4] = trial_mean

                erp = beta_mat @ v  # [times]
                erp_n.append(erp)
                meta.append(("conflict" if is_con else "no_conflict", tval))

            erp_n = np.column_stack(erp_n)
            n_rows = erp_n.shape[0] * erp_n.shape[1]
            T_time = np.tile(times, erp_n.shape[1])
            T_event = np.repeat("noun", n_rows)
            T_cond = np.repeat([m[0] for m in meta], erp_n.shape[0])
            T_mtf = np.repeat([m[1] for m in meta], erp_n.shape[0])
            T_chan = np.repeat(ch_name, n_rows)
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
        fix_times = np.unique(np.round(events.loc[events["type"] == "fixation", "trial_time"].dropna(), 1))
        fix_labels = ["elsewhere", "other", "target"]
        con_vals = [True, False]

        for ch_idx, ch_name in enumerate(chan_names):
            beta_mat = beta[ch_idx, :, :]
            erp_list = []
            meta = []
            for ft in fix_times:
                for cv in con_vals:
                    for fl in fix_labels:
                        # Construct covariate vector aligned to Unfold coefficient order.
                        # Only populate the fixation block (immediately after prev/next).
                        v = np.zeros(n_params, dtype=float)
                        fix_base = 2  # prev (1) + next (1) => fixation starts at index 2 (0-based)

                        conflict = 1.0 if cv else 0.0
                        is_other = 1.0 if fl == "other" else 0.0
                        is_target = 1.0 if fl == "target" else 0.0
                        ft = float(ft)

                        # fixation formula in Step G Julia script:
                        # 0 ~ 1 + condition * fix_at * trial_time + spl(saccAmpl, 5)
                        v[fix_base + 0] = 1.0
                        v[fix_base + 1] = conflict
                        v[fix_base + 2] = is_other
                        v[fix_base + 3] = is_target
                        v[fix_base + 4] = ft
                        v[fix_base + 5] = conflict * is_other
                        v[fix_base + 6] = conflict * is_target
                        v[fix_base + 7] = conflict * ft
                        v[fix_base + 8] = is_other * ft
                        v[fix_base + 9] = is_target * ft
                        v[fix_base + 10] = conflict * (is_other * ft)
                        v[fix_base + 11] = conflict * (is_target * ft)

                        # Spline terms: mirror the legacy MATLAB reconstruction behavior by
                        # setting all spline covariates to the mean saccade amplitude.
                        spline_start = fix_base + 12
                        if k_spline > 0:
                            v[spline_start:spline_start + k_spline] = mean_sacc

                        erp = beta_mat @ v  # [times]
                        erp_list.append(erp)
                        meta.append(("conflict" if cv else "no_conflict", fl, ft))

            if not erp_list:
                continue
            erp_f = np.column_stack(erp_list)
            n_rows = erp_f.shape[0] * erp_f.shape[1]
            T_time = np.tile(times, erp_f.shape[1])
            T_event = np.repeat("fixation", n_rows)
            T_cond = np.repeat([m[0] for m in meta], erp_f.shape[0])
            T_fixat = np.repeat([m[1] for m in meta], erp_f.shape[0])
            T_trtime = np.repeat([m[2] for m in meta], erp_f.shape[0])
            T_chan = np.repeat(ch_name, n_rows)
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
