from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
import warnings
import pandas as pd
from types import SimpleNamespace
import os

# Handles all data entering the simulation (command line arguments, taking data from the config file etc.). Doesn't encompass/ interact with conf.py though


def get_file_dir():
    return os.path.join(os.path.dirname(__file__))


def get_sim_configuration(config_name=None):
    df = pd.read_csv(
        get_file_dir() + "/configurations.csv", usecols=lambda x: x != "description"
    )
    df = df.set_index("name")
    print(df.index.tolist())
    if config_name is None:
        config_name = input("Choose the config: ")

    if config_name == "":
        config_name = "Default"
    return get_config_namespace_from_df_row(df, config_name)


def get_config_namespace_from_df_row(df: pd.DataFrame, config_name):
    config = df.loc[config_name]
    config = config.apply(pd.to_numeric)
    config["name"] = config.name
    return SimpleNamespace(**config)


def set_up_and_get_arguments():
    warnings.simplefilter("ignore")

    parser = ArgumentParser(formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument(
        "-sf", "--save_figure", help="Whether to save a graph of the simulation."
    )
    parser.add_argument(
        "-ow",
        "--overwrite",
        help="Whether to overwrite a pre-existing simulation if it exists.",
    )
    parser.add_argument(
        "-c",
        "--config",
        help="Name of the configuration in configurations.csv to load.",
    )
    args = vars(parser.parse_args())

    save_fig = args["save_figure"]
    overwrite = args["overwrite"]
    config_name = args["config"]

    return save_fig, overwrite, config_name
