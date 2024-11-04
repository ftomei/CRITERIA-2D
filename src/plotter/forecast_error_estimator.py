import os

import datetime
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.dates import DateFormatter
import matplotlib.dates as mdates
from sklearn.metrics import mean_squared_error, r2_score
import warnings

warnings.filterwarnings("ignore")



def main():
    obs_path = lambda x: os.path.join("/", "home", "data", "errano_all", "meteo" if x != "precipitation" else "water")
    forecast_path = os.path.join("/", "home", "data", "ground_truth_meteo")
    meteo_vars = ["air_humidity", "air_temperature", "solar_radiation", "wind_speed", "precipitation"]

    dfs = {
        meteo_var: pd.read_csv(os.path.join(obs_path(meteo_var), f"{meteo_var}.csv"))
        .set_index("timestamp")
        .join(
            other=pd.read_csv(os.path.join(forecast_path, f"{meteo_var}.csv")).set_index("timestamp"),
            lsuffix="_obs",
            rsuffix="_forecast",
            how="left"
            )
        # .apply(lambda row: 0 if row[f"{meteo_var}_obs"] == row[f"{meteo_var}_forecast"] else (abs(row[f"{meteo_var}_obs"] - row[f"{meteo_var}_forecast"]) / row[f"{meteo_var}_obs"]) * 100,
        #        axis=1)
        # .describe().loc[["mean", "std"]]
        for meteo_var in meteo_vars
    }

    for meteo_var in meteo_vars:
        error = mean_squared_error(
            dfs[meteo_var][f"{meteo_var}_forecast"],
            dfs[meteo_var][f"{meteo_var}_obs"],
            squared=False)
        print(f"{meteo_var}: {round(error, 2)}")


    # for meteo_var in meteo_vars:
    #     for var_type, value_path in {"obs": obs_path, "forecast": forecast_path}.items():
    #         for meteo_var in meteo_vars:
    #             load_path = os.path.join(value_path if var_type != "obs" else value_path(meteo_var), f"{meteo_var}.csv")
    #             dfs[var_type][meteo_var] = pd.read_csv(load_path)
    #             print(dfs[var_type][meteo_var])

    # for meteo_var in meteo_vars:
    #     for var_type in ["obs", "forecast"]:
    #         dfs[var_type][meteo_var] = pd.read_csv(load_path)



if __name__ == "__main__":
    main()