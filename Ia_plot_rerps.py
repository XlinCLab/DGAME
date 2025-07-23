import argparse
import time

from load_experiment import load_config, get_experiment_id, log_step_duration

def main(config: str | dict) -> dict:
    start_time = time.time()
    # Load experiment config
    if isinstance(config, str):
        config = load_config(config)
    experiment_id = get_experiment_id(config)


    # Log duration of this step in run config
    log_step_duration(config, start_time, step_id="Da_gaze_stats")

    return config


if __name__ == "__main__":
    parser = argparse.ArgumentParser("Calculate and plot gaze statistics.")
    parser.add_argument('config', help='Path to config.yml file')
    args = parser.parse_args()
    main(args.config)
