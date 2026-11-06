"""
Typed model of a DGAME experiment config, independent of any GUI toolkit.

The editable fields (and their types) are derived from a default config: every leaf value of the
default config becomes a `Field` whose type is taken from the default value's type.
"""

import copy
import os
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import yaml

from experiment import PARAM_ENABLED_KEY, RUN_CONFIG_KEY
from utils.run_config import load_config

CONFIG_DIR = Path(__file__).resolve().parents[2] / "config"
DEFAULTS_FILE_REGEX = re.compile(r"^dgame(\d+)_defaults\.yml$")
VERSION_KEY = "experiment.dgame_version"

# Top-level config sections which are not user-editable
HIDDEN_SECTIONS = (RUN_CONFIG_KEY, "config")

# Field kinds
BOOL, INT, FLOAT, STR, OPTIONAL_STR, STR_LIST, YAML = "bool", "int", "float", "str", "optional_str", "str_list", "yaml"
DIR_KEYS = {"root", "outdir"}


@dataclass
class Field:
    path: tuple[str, ...]
    default: Any
    kind: str
    help: str = ""

    @property
    def key(self) -> str:
        return ".".join(self.path)

    @property
    def label(self) -> str:
        return self.path[-1]

    @property
    def is_dir(self) -> bool:
        """Whether the field holds a directory path which can be chosen via a folder browser."""
        return self.label in DIR_KEYS

    @property
    def is_step_switch(self) -> bool:
        return len(self.path) == 4 and self.path[:2] == ("analysis", "steps") and self.label == PARAM_ENABLED_KEY


def infer_kind(value: Any) -> str:
    if isinstance(value, bool):
        return BOOL
    if isinstance(value, int):
        return INT
    if isinstance(value, float):
        return FLOAT
    if isinstance(value, str):
        return STR
    if value is None:
        return OPTIONAL_STR
    if isinstance(value, list) and all(isinstance(v, str) for v in value):
        return STR_LIST
    return YAML


def value_to_text(kind: str, value: Any) -> str:
    """Text representation of a value for display in a text-based editor."""
    if kind == STR_LIST:
        return "\n".join(str(v) for v in value)
    if kind == YAML:
        return yaml.safe_dump(value, default_flow_style=None, sort_keys=False).strip()
    return "" if value is None else str(value)


def text_to_value(kind: str, text: str, default: Any = None) -> Any:
    """Parse text from an editor into a value of the field's type. Raises ValueError if invalid."""
    if kind == INT:
        try:
            return int(text.strip())
        except ValueError:
            raise ValueError("expected a whole number") from None
    if kind == FLOAT:
        try:
            return float(text.strip())
        except ValueError:
            raise ValueError("expected a number") from None
    if kind == STR:
        return text
    if kind == OPTIONAL_STR:
        return text.strip() or None
    if kind == STR_LIST:
        return [line.strip() for line in text.splitlines() if line.strip()]
    if kind == YAML:
        try:
            value = yaml.safe_load(text) if text.strip() else None
        except yaml.YAMLError:
            raise ValueError("invalid YAML") from None
        if default is not None and not isinstance(value, type(default)):
            raise ValueError(f"expected a {type(default).__name__}")
        return value
    raise ValueError(f"cannot parse text for field of kind {kind}")


def extract_comments(yaml_text: str) -> dict[str, str]:
    """Map dotted config keys to their descriptions, taken from the comment lines directly above
    a key and from a trailing comment on its line."""
    comments: dict[str, str] = {}
    stack: list[tuple[int, str]] = []
    pending: list[str] = []
    for line in yaml_text.splitlines():
        stripped = line.strip()
        if not stripped:
            pending = []
            continue
        if stripped.startswith("#"):
            pending.append(stripped.lstrip("#").strip())
            continue
        match = re.match(r"^(\s*)([^\s:#\-][^:]*?|\"[^\"]*\"|'[^']*'):(?:\s+(.*))?$", line)
        if not match:
            pending = []
            continue
        indent, key, rest = len(match.group(1)), match.group(2).strip("\"'"), match.group(3) or ""
        while stack and stack[-1][0] >= indent:
            stack.pop()
        stack.append((indent, key))
        inline = re.search(r"\s#\s*(.*)$", rest)
        description = " ".join(pending + ([inline.group(1)] if inline else []))
        if description:
            comments[".".join(k for _, k in stack)] = description
        pending = []
    return comments


def defaults_file(version: str | int) -> Path:
    return CONFIG_DIR / f"dgame{version}_defaults.yml"


def available_versions() -> list[str]:
    """DGAME versions for which a default config file exists."""
    versions = [m.group(1) for p in CONFIG_DIR.iterdir() if (m := DEFAULTS_FILE_REGEX.match(p.name))]
    return sorted(versions, key=int)


def _comments_of_file(path: Path) -> dict[str, str]:
    """Comments of a config file, including those of any file it extends (`config.extends`)."""
    text = path.read_text()
    comments = {}
    extends = (yaml.safe_load(text) or {}).get("config", {}).get("extends")
    if extends is not None:
        base = Path(extends) if os.path.isabs(extends) else path.parent / extends
        comments.update(_comments_of_file(base))
    comments.update(extract_comments(text))
    return comments


def load_defaults(version: str | int) -> tuple[dict, dict[str, str]]:
    """Load the default config of a DGAME version and the descriptions of its parameters."""
    path = defaults_file(version)
    config = load_config(str(path))
    return config, _comments_of_file(path)


def flatten(config: dict, path: tuple[str, ...] = ()) -> list[tuple[tuple[str, ...], Any]]:
    """List the (path, value) leaves of a nested dict. Empty dicts count as leaves."""
    leaves = []
    for key, value in config.items():
        if isinstance(value, dict) and value:
            leaves.extend(flatten(value, path + (str(key),)))
        else:
            leaves.append((path + (str(key),), value))
    return leaves


def _get(config: dict, path: tuple[str, ...]) -> tuple[bool, Any]:
    node = config
    for key in path:
        if not isinstance(node, dict) or key not in node:
            return False, None
        node = node[key]
    return True, node


def _set(config: dict, path: tuple[str, ...], value: Any) -> None:
    node = config
    for key in path[:-1]:
        node = node.setdefault(key, {})
    node[path[-1]] = value


@dataclass
class ConfigModel:
    """Editable values of a DGAME config, typed after the default config they derive from."""
    defaults: dict
    comments: dict[str, str] = field(default_factory=dict)
    fields: dict[str, Field] = field(init=False)
    values: dict[str, Any] = field(init=False)
    errors: dict[str, str] = field(init=False)
    extras: dict = field(init=False)

    def __post_init__(self):
        visible = {k: v for k, v in self.defaults.items() if k not in HIDDEN_SECTIONS}
        self.defaults = copy.deepcopy(visible)
        self.fields = {}
        for path, default in flatten(self.defaults):
            f = Field(path=path, default=default, kind=infer_kind(default), help=self.comments.get(".".join(path), ""))
            self.fields[f.key] = f
        self.values = {k: copy.deepcopy(f.default) for k, f in self.fields.items()}
        self.errors = {}
        self.extras = {}

    # --- Loading and exporting ---

    def load(self, config: dict) -> None:
        """Set values from a (possibly partial) config. Parameters unknown to the defaults are preserved as extras."""
        self.extras = {}
        for path, value in flatten({k: v for k, v in config.items() if k not in HIDDEN_SECTIONS}):
            key = ".".join(path)
            if key in self.fields:
                self.values[key] = value
            elif any(key.startswith(k + ".") for k, f in self.fields.items() if f.kind == YAML):
                continue  # nested within a YAML-valued field; covered by that field's value
            else:
                _set(self.extras, path, value)
        # Values loaded into YAML fields given as non-empty dicts are flattened above; re-read them whole
        for key, f in self.fields.items():
            if f.kind == YAML:
                found, value = _get(config, f.path)
                if found:
                    self.values[key] = value

    def to_dict(self, only_changes: bool = False) -> dict:
        """Nested config dict of the current values (optionally only those differing from the defaults)."""
        result: dict = copy.deepcopy(self.extras)
        for key, f in self.fields.items():
            value = self.values[key]
            if only_changes and value == f.default and key != VERSION_KEY:
                continue
            _set(result, f.path, copy.deepcopy(value))
        return result

    def rebase(self, defaults: dict, comments: dict[str, str]) -> None:
        """Switch to new defaults (e.g. of another DGAME version), keeping the user's own edits."""
        edits = {k: v for k, v in self.values.items() if k in self.fields and v != self.fields[k].default}
        extras = self.extras
        self.defaults, self.comments = defaults, comments
        self.__post_init__()
        for key, value in edits.items():
            if key in self.fields and key != VERSION_KEY:
                self.values[key] = value
        self.extras = extras

    # --- Editing ---

    def set_text(self, key: str, text: str) -> str | None:
        """Set a field from editor text. Returns an error message, or None if the text was valid."""
        f = self.fields[key]
        try:
            self.values[key] = text_to_value(f.kind, text, f.default)
        except ValueError as exc:
            self.errors[key] = str(exc)
            return str(exc)
        self.errors.pop(key, None)
        return None

    def set_value(self, key: str, value: Any) -> None:
        self.values[key] = value
        self.errors.pop(key, None)

    # --- Validation ---

    def missing_required(self, required_patterns: list[str] | None) -> list[str]:
        """Keys of required fields which are empty."""
        missing = []
        for key, f in self.fields.items():
            if required_patterns and any(re.search(p, key) for p in required_patterns):
                value = self.values[key]
                if value is None or (isinstance(value, (str, list, dict)) and len(value) == 0):
                    missing.append(key)
        return missing

    def problems(self, required_patterns: list[str] | None = None) -> list[str]:
        """Human-readable problems which should be resolved before running the analysis."""
        problems = [f"{key}: {message}" for key, message in self.errors.items()]
        problems += [f"{key}: required" for key in self.missing_required(required_patterns)]
        root = self.values.get("experiment.input_data.root")
        if root and not os.path.isdir(root):
            problems.append("experiment.input_data.root: folder does not exist")
        return problems


def discover_subjects(model: ConfigModel) -> list[str]:
    """Subject IDs found as subdirectories of the input data's recordings/xdf directory."""
    values = model.values
    root = values.get("experiment.input_data.root")
    if not root:
        return []
    xdf_dir = os.path.join(root, values.get("experiment.input_data.recordings_dir", ""), values.get("experiment.input_data.xdf_dir", ""))
    if not os.path.isdir(xdf_dir):
        return []
    return sorted(d for d in os.listdir(xdf_dir) if os.path.isdir(os.path.join(xdf_dir, d)) and not d.startswith("."))


def dependency_status() -> dict[str, bool]:
    """Whether the external software which some analysis steps depend upon was found."""
    import shutil
    return {"R": shutil.which("Rscript") is not None, "Julia": shutil.which("julia") is not None}
