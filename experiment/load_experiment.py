import json
import logging
import os
import time
from datetime import timedelta
from typing import Self, Union

from experiment.constants import RUN_CONFIG_KEY
from experiment.test_subjects import parse_subject_ids
from utils.run_config import dump_config, load_config
from utils.utils import create_timestamp

logger = logging.getLogger(__name__)


class Experiment:
    def __init__(self,
                 config_path: str,
                 ):
        self.start_time = time.time()
        self.config = self.load_config(config_path)
        self.experiment_id = self.get_experiment_id()
        self.outdir = self.create_experiment_outdir()
        self.subjects = self.get_parameter("subjects")
        self.subject_ids, self.subject_id_regex = parse_subject_ids(self.subjects)

    def load_config(self, config: str | dict):
        """Load experiment's run config."""
        if isinstance(config, str):
            config = os.path.abspath(config)
            logger.info(f"Initializing experiment from {config} ...")
            config = load_config(config)
            logger.info(json.dumps(config, indent=4))
            return config
        elif isinstance(config, dict):
            return config
        else:
            raise TypeError(f"Expected config to be a filepath string or dict, found {type(config)}")

    @classmethod
    def from_input(cls: type[Self],
                   experiment: Union[Self, str, dict]) -> Self:
        """Initialize input as an Experiment (sub)class object."""
        if isinstance(experiment, cls):
            return experiment
        return cls(experiment)

    def get_parameter(self, *parameter_keys: str, default = None):
        """Retrieve a parameter value from the experiment config."""
        value = self.config["experiment"]
        try:
            for key in parameter_keys:
                value = value[key]
            return value
        except (KeyError, TypeError):
            return default

    def get_experiment_id(self, add_timestamp: bool = False) -> str:
        """Retrieve experiment ID from config (if set) and optionally combine with timestamp."""
        experiment_id = self.get_parameter("id")
        _, timestamp = create_timestamp()
        if experiment_id is None or experiment_id.strip() == "":
            experiment_id = timestamp
        elif add_timestamp:
            experiment_id = os.path.join(experiment_id, timestamp)
        self.config[RUN_CONFIG_KEY]["id"] = experiment_id
        return experiment_id

    def create_experiment_outdir(self) -> str:
        """Retrieve and create experiment output directory."""
        base_output_dir = self.get_parameter("outdir")
        if base_output_dir is None or base_output_dir.strip() == "":
            base_output_dir = "out"
        experiment_id = self.experiment_id
        output_dir = os.path.join(os.path.abspath(base_output_dir), experiment_id)
        if os.path.exists(output_dir):
            # If experiment outdir already exists, ask user to confirm
            # in order to avoid potentially overwriting previous results
            def ask_user_to_confirm_overwrite(output_dir):
                valid_answers = {"y", "yes", "n", "no"}
                overwrite_previous = None
                while overwrite_previous not in valid_answers:
                    overwrite_previous = input(f"Output directory {output_dir} already exists. Overwrite previous results? [Y/N]")
                    overwrite_previous = overwrite_previous.lower().strip()
                return overwrite_previous in {"y", "yes"}
            overwrite_previous = ask_user_to_confirm_overwrite(output_dir)
            if overwrite_previous is False:
                # If user specifies to NOT overwrite previous results, add timestamp
                # in order to disambiguate and keep previous results
                _, timestamp = create_timestamp()
                experiment_id = "_".join([experiment_id, timestamp])
                output_dir = os.path.join(os.path.abspath(base_output_dir), experiment_id)

        os.makedirs(output_dir, exist_ok=True)
        logger.info(f"Created experiment output directory: {output_dir}")
        self.config[RUN_CONFIG_KEY]["outdir"] = output_dir
        return output_dir

    def log_step_duration(self, start_time: float, step_id: str) -> str:
        """Calculate and log duration of a particular experiment processing step."""
        end_time = time.time()
        duration = str(timedelta(seconds=int(end_time - start_time)))
        self.config[RUN_CONFIG_KEY]["duration"][step_id] = duration
        logger.info(f"Step {step_id} completed successfully (duration: {duration}).")
        return duration

    def log_total_duration(self):
        """Calculate total experiment duration."""
        end_time = time.time()
        total_duration = str(timedelta(seconds=int(end_time - self.start_time)))
        self.config[RUN_CONFIG_KEY]["duration"]["total"] = total_duration
        logger.info(f"Total experiment duration: {total_duration}")
    
    def finish(self):
        """Log total duration and write updated experiment config to output directory."""
        self.log_total_duration()
        config_outpath = os.path.join(self.outdir, "config.yml")
        dump_config(self.config, config_outpath)
        logger.info(f"Wrote experiment run config to {os.path.abspath(config_outpath)}")
