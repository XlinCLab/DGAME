import pandas as pd
import rpy2.robjects as robjects
from rpy2.robjects import StrVector, pandas2ri
from rpy2.robjects import r as r_interface
from rpy2.robjects.vectors import DataFrame as RDataFrame

from utils import generate_variable_name

# Source R functions from dependencies.R and load and/or install dependencies
robjects.r["source"]("dependencies.R")


def r_install_package(package: str) -> None:
    """Installs a single R package."""
    robjects.globalenv["install_if_needed"](package)


def r_install_packages(package_list: list) -> None:
    """Installs multiple R packages."""
    robjects.globalenv["install_packages_if_needed"](StrVector(package_list))


def r_eval(expression: str, name: str = None):
    """Evaluate an expression in R and assign it to a variable, then retrieve and return that variable from the R environment."""
    if name is None:
        # Generate random variable name string if none provided
        name = generate_variable_name()
    r_interface(f"{name} <- {expression}")
    result = r_interface[name][0]
    return result


def convert_pandas2r_dataframe(pd_dataframe: pd.DataFrame) -> RDataFrame:
    """Converts a pandas dataframe to an R dataframe."""
    with (robjects.default_converter + pandas2ri.converter).context():
        r_dataframe = pandas2ri.py2rpy(pd_dataframe)
        return r_dataframe


def convert_r2pandas_dataframe(r_dataframe: RDataFrame) -> pd.DataFrame:
    """Converts an R dataframe to a pandas dataframe."""
    with (robjects.default_converter + pandas2ri.converter).context():
        pd_dataframe = pandas2ri.rpy2py(r_dataframe)
        return pd_dataframe
