import os
import yaml

def recursively_inherit_dict_values(target: dict, source: dict) -> None:
    for key, value in source.items():
        if target is not None and key not in target:
            target[key] = value
        elif isinstance(value, dict):
            recursively_inherit_dict_values(target[key], value)

def to_absolute_path(path: str) -> str:
    if os.path.isabs(path):
        return path
    current_file_path = os.path.abspath(__file__)
    root_path = os.path.abspath(os.path.join(current_file_path, os.pardir))
    return os.path.join(root_path, path)


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

    included_config_block: dict | None = loaded_config.get('config')
    if included_config_block is None:
        return loaded_config
    included_config_path: str | None = included_config_block.get('extend')
    if included_config_path is None:
        return loaded_config

    included_config = None
    if included_config_path is not None:
        included_config_path: str = to_absolute_path(included_config_path)
        included_config = load_config(included_config_path)

    recursively_inherit_dict_values(loaded_config, included_config)

    return loaded_config