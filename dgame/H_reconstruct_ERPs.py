import argparse
import os
from typing import Dict

import numpy as np
import pandas as pd

from dgame.constants import STEP_H_KEY
from experiment.load_experiment import Experiment


def _load_events(events_file: str) -> pd.DataFrame:
    df = pd.read_csv(events_file)
    if "trial_time" in df.columns:
        df["trial_time"] = pd.to_numeric(df["trial_time"], errors="coerce")
    if "trial" in df.columns:
        df["trial"] = pd.to_numeric(df["trial"], errors="coerce")
    if "saccAmpl" in df.columns:
        df["saccAmpl"] = pd.to_numeric(df["saccAmpl"], errors="coerce")
    return df


def _make_spline_basis(x: float, df: int = 5) -> np.ndarray:
    try:
        import patsy

        mat = patsy.dmatrix("bs(x, df=%d, include_intercept=False)" % df, {"x": [x]}, return_type="dataframe")
        return np.asarray(mat, dtype=float).reshape(-1)
    except Exception:
        return np.array([x ** (i + 1) for i in range(df)], dtype=float)


def _vector_from_params(param_names: list[str], values: Dict[str, float]) -> np.ndarray:
    v = np.zeros(len(param_names), dtype=float)
    for name, val in values.items():
        if name in param_names:
            v[param_names.index(name)] = val
    return v


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

        uf_file = os.path.join(unfold_out_dir, f"{subject_id}_ufresult.npz")
        events_file = os.path.join(subject_eeg_dir, f"{subject_id}_director_events.csv")

        uf = np.load(uf_file, allow_pickle=True)
        beta = uf["beta"]  # shape (n_channels, n_times, n_params)
        times = uf["times"] * 1000.0
        param_names = [str(x) for x in uf["param_names"]]
        chan_names = [str(x) for x in uf["chan_names"]]

        events = _load_events(events_file)
        mean_tg = np.nanmean(events.loc[events["type"] == "N", "trial_time"].to_numpy())
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
            beta_mat = beta[ch_idx, :, :]
            if all_times_n.size == 0:
                continue

            erp_n = []
            meta = []
            for tval, is_con in zip(all_times_n, is_con_n):
                if is_con:
                    values = {
                        "N:Intercept": 1.0,
                        "N:condition_conflict": 1.0,
                        "N:trial": trial_con,
                        "N:trial_time": tval,
                        "N:condition_conflict:trial_time": tval,
                    }
                    cond_label = "conflict"
                else:
                    values = {
                        "N:Intercept": 1.0,
                        "N:condition_conflict": 0.0,
                        "N:trial": trial_nocon,
                        "N:trial_time": tval,
                        "N:condition_conflict:trial_time": 0.0,
                    }
                    cond_label = "no_conflict"
                v = _vector_from_params(param_names, values)
                erp = beta_mat @ v
                erp_n.append(erp)
                meta.append((cond_label, tval))

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
                        values = {
                            "fixation:Intercept": 1.0,
                            "fixation:condition_conflict": 1.0 if cv else 0.0,
                            "fixation:fix_at_other": 1.0 if fl == "other" else 0.0,
                            "fixation:fix_at_target": 1.0 if fl == "target" else 0.0,
                            "fixation:trial_time": ft,
                            "fixation:condition_conflict:fix_at_other": (1.0 if cv and fl == "other" else 0.0),
                            "fixation:condition_conflict:fix_at_target": (1.0 if cv and fl == "target" else 0.0),
                            "fixation:condition_conflict:trial_time": (ft if cv else 0.0),
                            "fixation:fix_at_other:trial_time": (ft if fl == "other" else 0.0),
                            "fixation:fix_at_target:trial_time": (ft if fl == "target" else 0.0),
                            "fixation:condition_conflict:fix_at_other:trial_time": (ft if cv and fl == "other" else 0.0),
                            "fixation:condition_conflict:fix_at_target:trial_time": (ft if cv and fl == "target" else 0.0),
                        }
                        # add spline terms if present
                        spline_vals = _make_spline_basis(mean_sacc, df=5)
                        for i in range(1, 6):
                            key = f"fixation:saccAmpl_spline_{i}"
                            if key in param_names and i - 1 < len(spline_vals):
                                values[key] = float(spline_vals[i - 1])
                        v = _vector_from_params(param_names, values)
                        erp = beta_mat @ v
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
