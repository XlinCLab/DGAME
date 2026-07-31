import os

import yaml

from experiment import RUN_CONFIG_KEY
from utils.utils import (convert_sets_to_lists, create_timestamp,
                         get_git_commit_hash, recursively_inherit_dict_values)


def populate_config_from_extended(config_path: str, visited: frozenset[str] = frozenset()) -> dict:
    """Load a YAML config file and recursively merge in its `config.extends` chain (if any),
    without attaching any run-specific metadata. Each file in the chain fills in only the 
    values left unspecified by the ones that extend it."""
    config_path = os.path.abspath(config_path)
    if config_path in visited:
        chain = " -> ".join(list(visited) + [config_path])
        raise ValueError(f"Circular config.extends chain detected: {chain}")

    with open(config_path, 'r') as f:
        loaded_config: dict = yaml.safe_load(f)

    extends = loaded_config.get("config", {}).get("extends")
    if extends is not None:
        if not os.path.isabs(extends):
            extends = os.path.abspath(os.path.join(os.path.dirname(config_path), extends))
        base_config = populate_config_from_extended(extends, visited=visited | {config_path})
        recursively_inherit_dict_values(loaded_config, base_config)

    return loaded_config


def load_config(config_path: str, default_config: str | dict = None) -> dict:
    """Returns a dictionary containing parameters from a specified config.yml file, recursively
    merged with any `config.extends` base file(s) it declares.

    Args:
        config_path (str): Path to config.yml file
        default_config (str | dict): Optional explicit default config (path or already-loaded
            dict) to additionally inherit unspecified values from. Takes precedence over (i.e.
            is applied on top of, filling in only what's still missing after) any config.extends
            chain resolved from config_path itself.

    Returns:
        config: nested dictionary of parameter names and values
    """
    loaded_config = populate_config_from_extended(config_path)

    # Recursively inherit unspecified values from default config
    if default_config is not None:
        if isinstance(default_config, str):
            default_config = populate_config_from_extended(default_config)
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
