import pandas as pd
import rpy2.robjects as robjects
from rpy2.robjects import StrVector, globalenv, pandas2ri

# Source R functions from dependencies.R and load and/or install dependencies
robjects.r["source"]("dependencies.R")


def r_install_package(package: str) -> None:
    """Installs a single R package."""
    robjects.globalenv["install_if_needed"](package)


def r_install_packages(package_list: list) -> None:
    """Installs multiple R packages."""
    robjects.globalenv["install_packages_if_needed"](StrVector(package_list))


def convert_pandas2r_dataframe(pd_dataframe: pd.DataFrame) -> robjects.vectors.DataFrame:
    """Converts a pandas dataframe to an R dataframe."""
    with (robjects.default_converter + pandas2ri.converter).context():
        r_dataframe = pandas2ri.py2rpy(pd_dataframe)
        return r_dataframe
