import glob
import logging
import os
import re
from collections import defaultdict
from typing import Iterable

logger = logging.getLogger(__name__)


def parse_subject_ids(subject_ids: Iterable | str | int | None) -> tuple[list, str]:
    """Parse subject IDs from string, iterable, or integer and return a list of parsed IDs and a regular expression to match them."""
    if len(subject_ids) == 0 or subject_ids is None or isinstance(subject_ids, str) and subject_ids.strip() == "":
        subject_id_regex = r"*"
        subject_ids = []
    elif isinstance(subject_ids, str):
        subject_id_regex = subject_ids.strip()
        subject_ids = [subject_ids]
    elif isinstance(subject_ids, int):
        subject_id_regex = str(subject_ids)
        subject_ids = [str(subject_ids)]
    else:
        subject_id_regex = "(" + "|".join([str(s) for s in subject_ids]) + ")"
    return subject_ids, subject_id_regex


def list_subject_files(dir: str,
                       subject_regex: str | None = None,
                       suffix: str | None = None,
                       recursive: bool = False) -> list:
    """Filter a directory and returns a list of files matching a subject ID and optional file suffix."""
    # Ensure that the base directory exists
    dir = os.path.abspath(dir)
    if not os.path.exists(dir):
        raise FileNotFoundError(f"Directory not found: {dir}")
    # Compile a regular expression matching subject(s) and relevant suffix
    suffix = r"*" if suffix is None else suffix
    subject_regex = "" if subject_regex is None else subject_regex
    path_regex = re.compile(f"{subject_regex}{suffix}$")
    # Filter files in directory by regex
    subject_files = [
        os.path.join(dir, filepath) for filepath in glob.glob("**", root_dir=dir, recursive=recursive)
        if path_regex.match(os.path.basename(filepath))
    ]
    if len(subject_files) == 0:
        logger.warning(f"No matching subject files found in {dir}")
    return sorted(subject_files)


def subject_files_dict(dir: str,
                       subject_regex: str | None = None,
                       suffix: str | None = None,
                       recursive: bool = False) -> defaultdict:
    """Return a dictionary of subject IDs and files per subject ID matching a given pattern within a directory."""
    subject_files = list_subject_files(
        dir=dir, subject_regex=subject_regex, suffix=suffix, recursive=recursive,
    )
    suffix = r"*" if suffix is None else suffix
    files_per_subject = defaultdict(lambda: [])
    for subject_file in subject_files:
        subject_id = re.search(rf"{subject_regex}(?={suffix})", os.path.basename(subject_file))
        if not subject_id:
            logger.error(f"Could not extract subject ID from {subject_file}")
            continue
        subject_id = subject_id.group()
        files_per_subject[subject_id].append(subject_file)
    return files_per_subject


def subject_dirs_dict(root_dir: str,
                      subject_regex: str | None = None,
                      ) -> dict:
    """List subdirectories matching a subject regex."""
    subject_regex = r"*" if subject_regex is None else subject_regex
    subject_regex = re.compile(subject_regex)
    subject_subdirs = defaultdict(lambda: [])
    for path in glob.glob(os.path.join(root_dir, "*")):
        if os.path.isdir(path):
            subject_id = subject_regex.match(os.path.basename(path))
            if subject_id:
                subject_subdirs[subject_id.group()].append(path)
    return subject_subdirs
