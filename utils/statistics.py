import numpy as np
import pandas as pd
from statsmodels.regression.linear_model import RegressionResultsWrapper
from statsmodels.stats.multitest import multipletests

ALPHA = 0.05  # p-value significance threshold for statistics tests


def summarize_stats_model(model: RegressionResultsWrapper) -> pd.DataFrame:
    """Return summary dataframe for a fitted regression model."""
    model_summary = pd.DataFrame({
        'predictor': model.params.index,
        'coef': model.params.values,
        'p_value': model.pvalues.values,
        't_value': model.tvalues.values,
    })
    return model_summary


def fdr_adjust_pvals(p_values: np.ndarray, alpha: float = ALPHA):
    """Run FDR correction (Benjamini-Hochberg) for an array of p-values from permutation testing."""
    _, pvals_corrected, _, _ = multipletests(
        p_values,
        alpha=alpha,
        method='fdr_bh',  # equivalent to R's method = "fdr" in p.adjust
    )
    return pvals_corrected
