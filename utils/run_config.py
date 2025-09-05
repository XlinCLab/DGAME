import yaml

from experiment.constants import RUN_CONFIG_KEY
from utils.utils import (convert_sets_to_lists, create_timestamp,
                         get_git_commit_hash, recursively_inherit_dict_values)


def load_config(config_path: str, default_config: str | dict = None) -> dict:
    """Returns a dictionary containing parameters from a specified config.yml file

    Args:
        config_path (str): Path to config.yml file

    Returns:
        config: nested dictionary of parameter names and values
    """
    with open(config_path, 'r') as f:
        loaded_config: dict = yaml.safe_load(f)

    # Recursively inherit unspecified values from default config
    if default_config is not None:
        if isinstance(default_config, str):
            default_config = load_config(default_config, default_config=None)
        else:
            assert isinstance(default_config, dict)
        recursively_inherit_dict_values(loaded_config, default_config)

    return init_run_config(loaded_config)


def init_run_config(config: dict):
    """Initialize run config with software version, start time, and duration."""
    _, timestamp = create_timestamp()
    run_config = {
        "version": get_git_commit_hash(),
        "start_time": timestamp,
        "duration": {},
    }
    config[RUN_CONFIG_KEY] = run_config

    return config


def dump_config(config, config_outpath):
    """Dumps a YAML config to an output file."""
    with open(config_outpath, "w") as f:
        yaml.dump(convert_sets_to_lists(config), f)
