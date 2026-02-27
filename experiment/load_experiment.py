import json
import logging
import os
import time
from collections import defaultdict
from datetime import timedelta
from typing import Callable, Self, Union

from experiment.constants import RUN_CONFIG_KEY
from experiment.test_subjects import (parse_subject_ids, subject_dirs_dict,
                                      subject_files_dict)
from utils.run_config import dump_config, load_config
from utils.utils import (LogMessageCollector, create_timestamp,
                         recursively_inherit_dict_values)


class Experiment:
    def __init__(self,
                 config_path: str,
                 default_config: str = None,
                 logger: logging.Logger = None,
                 log_file: str = None,
                 ):
        self.start_time = time.time()
        self.logging_handler = LogMessageCollector()
        self.logger = self.configure_logger(logger)
        self.defaults = self.load_config(default_config) if default_config is not None else default_config
        self.config = self.load_config(config_path, default_config=self.defaults)
        self.experiment_id = self.get_experiment_id()
        self.outdir = self.create_experiment_outdir()
        self.logdir = self.create_experiment_logs_dir()
        self.experiment_logfile, self.logfile_handler = self.set_experiment_logfile(log_file=log_file)
        self.subjects = self.get_experiment_parameter("subjects")
        self.subject_ids, self.subject_id_regex = parse_subject_ids(self.subjects)

    def configure_logger(self,
                         logger: logging.Logger = None,
                         ) -> logging.Logger:
        """Configure the experiment logger."""
        if logger is None:
            logger = logging.getLogger(__name__)
        
        # Add custom logging handler
        logger.addHandler(self.logging_handler)
        return logger

    def load_config(self, config: str | dict, default_config: str | dict = None):
        """Load experiment's run config."""
        if isinstance(config, str):
            config = os.path.abspath(config)
            self.logger.info(f"Initializing experiment from {config} ...")
            if default_config and isinstance(default_config, str):
                default_config = load_config(config)
            config = load_config(config, default_config=default_config)
            self.logger.info(json.dumps(config, indent=4))
            return config
        elif isinstance(config, dict):
            if default_config:
                recursively_inherit_dict_values(config, default_config)
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

    def get_parameter(self, *parameter_keys: str, default=None):
        """Retrieve a parameter value from the run config."""
        value = self.config
        try:
            for key in parameter_keys:
                value = value[key]
            return value
        except (KeyError, TypeError):
            return default

    def get_experiment_parameter(self, *parameter_keys: str, default=None):
        """Retrieve a parameter value from the experiment config."""
        return self.get_parameter("experiment", *parameter_keys, default=default)

    def get_input_data_path(self, *parameter_keys: str, default=None):
        """Retrieve a data path from the experiment's input data config."""
        return self.get_parameter("experiment", "input_data", *parameter_keys, default=default)

    def get_analysis_parameter(self, *parameter_keys: str, default=None):
        """Retrieve a parameter value from the analysis config."""
        return self.get_parameter("analysis", *parameter_keys, default=default)

    def get_experiment_id(self, add_timestamp: bool = False) -> str:
        """Retrieve experiment ID from config (if set) and optionally combine with timestamp."""
        experiment_id = self.get_experiment_parameter("id")
        _, timestamp = create_timestamp()
        if experiment_id is None or experiment_id.strip() == "":
            experiment_id = timestamp
        elif add_timestamp:
            experiment_id = os.path.join(experiment_id, timestamp)
        self.config[RUN_CONFIG_KEY]["id"] = experiment_id
        return experiment_id

    def create_experiment_outdir(self) -> str:
        """Retrieve and create experiment output directory."""
        base_output_dir = self.get_experiment_parameter("outdir")
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
        self.logger.info(f"Created experiment output directory: {output_dir}")
        self.config[RUN_CONFIG_KEY]["outdir"] = output_dir
        return output_dir

    def create_experiment_logs_dir(self) -> str:
        """Create experiment logs directory."""
        logs_directory = os.path.join(self.outdir, "logs")
        os.makedirs(logs_directory, exist_ok=True)
        return logs_directory

    def set_experiment_logfile(self,
                               log_file: str = "experiment.log",
                               ) -> tuple[str, logging.FileHandler]:
        """
        Configure logging to experiment log file.
        Returns a tuple of the string path to the log file and its FileHandler.
        """
        log_file = os.path.abspath(os.path.join(self.logdir, log_file))
        self.logger.info(f"Experiment run log file: {log_file}")

        # Add any preceding log messages (from early experiment initialization) to this log file
        log_messages = self.logging_handler.messages
        with open(log_file, "w") as f:
            f.write("\n".join(log_messages))

        fh = logging.FileHandler(log_file)
        fh.setLevel(logging.INFO)
        formatter = logging.Formatter("%(asctime)s %(levelname)s: %(message)s")
        fh.setFormatter(formatter)
        self.logger.addHandler(fh)

        return log_file, fh

    def get_subject_dirs_dict(self, root_dir: str) -> dict:
        """Return a dictionary of per-subject directories from a root directory."""
        return subject_dirs_dict(
            root_dir=root_dir,
            subject_regex=self.subject_id_regex,
        )

    def get_subject_files_dict(self,
                               dir: str,
                               suffix: str | None = None,
                               recursive: bool = False
                               ) -> defaultdict:
        """Return a dictionary of per-subject files within a directory."""
        return subject_files_dict(
            dir=dir,
            subject_regex=self.subject_id_regex,
            suffix=suffix,
            recursive=recursive,
        )

    def log_step_duration(self, start_time: float, step_id: str) -> str:
        """Calculate and log duration of a particular experiment processing step."""
        end_time = time.time()
        duration = str(timedelta(seconds=int(end_time - start_time)))
        self.config[RUN_CONFIG_KEY]["duration"][step_id] = duration
        self.logger.info(f"Step {step_id} completed successfully (duration: {duration}).")
        return duration

    def log_total_duration(self):
        """Calculate total experiment duration."""
        end_time = time.time()
        total_duration = str(timedelta(seconds=int(end_time - self.start_time)))
        self.config[RUN_CONFIG_KEY]["duration"]["total"] = total_duration
        self.logger.info(f"Total experiment duration: {total_duration}")

    def finish(self):
        """Log total duration and write updated experiment config to output directory."""
        self.log_total_duration()
        config_outpath = os.path.join(self.outdir, "config.yml")
        dump_config(self.config, config_outpath)
        self.logger.info(f"Wrote experiment run config to {os.path.abspath(config_outpath)}")


class ExperimentStep:
    def __init__(self,
                 label: str,
                 main_func: Callable,
                 experiment: Experiment,
                 log_file: str = None,
                 ):
        self.label = label
        self.experiment = experiment
        self.logfile_handler = None
        self.log_file = self.set_logfile(log_file)
        self.logger = self.configure_logger(experiment.logger)
        self.main_func = main_func
        self.start_time = None
        self.duration = None
    
    def set_logfile(self, log_file: str = None) -> str:
        if log_file is None:
            log_file = os.path.join(self.experiment.logdir, f"{self.label}.log")
        log_dir = os.path.dirname(log_file)
        os.makedirs(log_dir, exist_ok=True)
        return log_file
    
    def configure_logger(self, logger: logging.Logger) -> logging.Logger:
        """Configure logging to log file."""
        fh = logging.FileHandler(self.log_file)
        fh.setLevel(logging.INFO)
        formatter = logging.Formatter("%(asctime)s %(levelname)s: %(message)s")
        fh.setFormatter(formatter)
        logger.addHandler(fh)
        self.logfile_handler = fh
        return logger

    def run(self):
        self.logger.info(f"Running experiment component {self.label} ...")
        self.logger.info(f"Log file: {self.log_file}")
        self.start_time = time.time()
        result = self.main_func(self.experiment)
        self.log_duration()
        self.logger.removeHandler(self.logfile_handler)
        return result

    def log_duration(self) -> str:
        self.duration = self.experiment.log_step_duration(
            self.start_time,
            step_id=self.label
        )
        return self.duration
