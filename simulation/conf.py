import pandas as pd
from types import SimpleNamespace


interrupt = False
debug = True
path_to_data = "..sim_data/sim.pickle"
path_to_output = "..outputs/"
m_adjustment = True  # Make my simulation seem more like Marta's (But this can likely be removed, since it's just for small temporary changes)
sim_state_init_type = "detailed"


def get_sim_configuration(config_name=None):
    df = pd.read_csv("conf.csv", usecols=lambda x: x != "description")
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
