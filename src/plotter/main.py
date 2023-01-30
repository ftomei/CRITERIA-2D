import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import os
import time


def main():

    tuning_folder = os.path.join("data", "errano_tuning")
    evaluation_folder = os.path.join("data", "errano_evaluation_aug")
    output_folder =  os.path.join("plots")

    meteo_dict = {}
    meteo_vars = {
        "air_humidity": "%",
        "air_temperature": "°C",
        "solar_radiation": "W/mq",
        "wind_speed": "m/s",
    }

    for meteo_var in meteo_vars.keys():
        meteo_dict[meteo_var] = (
            pd.read_csv(os.path.join(tuning_folder, "meteo", f"{meteo_var}.csv"))
            .append(
                pd.read_csv(
                    os.path.join(evaluation_folder, "meteo", f"{meteo_var}.csv")
                ),
                ignore_index=True,
            )
            .set_index("timestamp")
        )

    meteo_df = pd.concat(meteo_dict.values(), axis=1)
    meteo_df.index = pd.to_datetime(meteo_df.index, unit="s")

    print(meteo_df)

    # define subplot layout
    nrows, ncols = 2, 2
    fig, axes = plt.subplots(nrows=nrows, ncols=ncols)

    # add DataFrames to subplots

    for idx, meteo_var in reversed(list(enumerate(meteo_vars.keys()))):
        meteo_df[meteo_var].plot(ax=axes[int(idx / nrows), idx % nrows], sharex=True)
        axes[int(idx / nrows), idx % nrows].set_xlim(
            [meteo_df.index[0], meteo_df.index[-1]]
        )
        axes[int(idx / nrows), idx % nrows].set_ylim(
            [0, 1000] if meteo_var == "solar_radiation" else (
                [0, 15] if meteo_var == "wind_speed" else (
                [0, 100] if meteo_var == "air_humidity" else [0, 40]))
        )
        axes[int(idx / nrows), idx % nrows].set_xlabel("")
        if meteo_var == "wind_speed":
            xticks = axes[int(idx / nrows), idx % nrows].get_xticks()
        axes[int(idx / nrows), idx % nrows].set_xticks(xticks)
        if idx < 2:
            axes[int(idx / nrows), idx % nrows].tick_params(length=0)
        axes[int(idx / nrows), idx % nrows].set_ylabel(meteo_vars[meteo_var], rotation=0, labelpad=13)
        axes[int(idx / nrows), idx % nrows].set_title(
            meteo_var.replace("_", " ").capitalize()
        )

    plt.tight_layout()
    fig.set_size_inches(13, 6)
    fig.savefig(os.path.join(output_folder, "meteo.pdf"))
    fig.savefig(os.path.join(output_folder, "meteo.png"))


if __name__ == "__main__":
    main()
