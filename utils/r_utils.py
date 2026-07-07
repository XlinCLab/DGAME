import pandas as pd

from utils.utils import generate_variable_name

# NB: rpy2 imports are deferred into each function rather than done at module level,
# since importing rpy2 embeds an R interpreter in-process and mutates LD_LIBRARY_PATH,
# which can break subprocess-based tools (e.g. Julia) if that happens too early


def r_eval(expression: str, name: str = None):
    """Evaluate an expression in R and assign it to a variable, then retrieve and return that variable from the R environment."""
    from rpy2.robjects import r as r_interface
    from rpy2.robjects.vectors import DataFrame as RDataFrame
    from rpy2.robjects.vectors import Vector

    if name is None:
        # Generate random variable name string if none provided
        name = generate_variable_name()
    r_interface(f"{name} <- {expression}")
    result = r_interface[name]

    # Return first element if scalar, full object otherwise
    if isinstance(result, Vector) and not isinstance(result, RDataFrame) and len(result) == 1:
        return result[0]
    return result


def convert_pandas2r_dataframe(pd_dataframe: pd.DataFrame):
    """Converts a pandas dataframe to an R dataframe."""
    import rpy2.robjects as robjects
    from rpy2.robjects import pandas2ri
    with (robjects.default_converter + pandas2ri.converter).context():
        r_dataframe = pandas2ri.py2rpy(pd_dataframe)
        return r_dataframe


def convert_r2pandas_dataframe(r_dataframe) -> pd.DataFrame:
    """Converts an R dataframe to a pandas dataframe."""
    import rpy2.robjects as robjects
    from rpy2.robjects import pandas2ri
    with (robjects.default_converter + pandas2ri.converter).context():
        pd_dataframe = pandas2ri.rpy2py(r_dataframe)
        return pd_dataframe
