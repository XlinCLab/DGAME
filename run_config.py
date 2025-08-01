import yaml

from constants import RUN_CONFIG_KEY
from utils import (convert_sets_to_lists, create_timestamp,
                   get_git_commit_hash, to_absolute_path)


def load_config(config_path: str) -> dict:
    """Returns a dictionary containing parameters from a specified config.yml file

    Args:
        config_path (str): Path to config.yml file

    Returns:
        config: nested dictionary of parameter names and values
    """

    config_path = to_absolute_path(config_path)
    with open(config_path, 'r') as f:
        loaded_config: dict = yaml.safe_load(f)

    # TODO add default config from which to recursively inherit values
    # recursively_inherit_dict_values(loaded_config, included_config)

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
