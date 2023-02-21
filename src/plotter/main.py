import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.dates import DateFormatter
import matplotlib.dates as mdates
from sklearn.metrics import mean_squared_error
import warnings
warnings.filterwarnings("ignore")

import os
import time

forbidden_sensors = ["z20_y0_x0", "z20_y0_x25", "z40_y0_x0", "z40_y0_x50"]

def meteo(
    tuning_folder=os.path.join("data", "errano_tuning"),
    evaluation_folder=os.path.join("data", "errano_evaluation_1gg"),
    output_folder=os.path.join("plots"),
):
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

    df = pd.concat(meteo_dict.values(), axis=1)
    df = df.reset_index()
    df["timestamp"] = pd.to_datetime(df["timestamp"], unit="s")
    # df["timestamp"] = df["timestamp"].dt.strftime('%Y-%m-%d')
    df = df.set_index("timestamp")

    # print(df)

    # define subplot layout
    nrows, ncols = 2, 2
    fig, axes = plt.subplots(nrows=nrows, ncols=ncols)

    # add DataFrames to subplots

    for idx, meteo_var in reversed(list(enumerate(meteo_vars.keys()))):
        df[meteo_var].plot(ax=axes[int(idx / nrows), idx % nrows], sharex=True)
        # axes[int(idx / nrows), idx % nrows].set_xlim(
        #     [0, df.shape[0]]
        # )
        axes[int(idx / nrows), idx % nrows].set_ylim(
            [0, 1000]
            if meteo_var == "solar_radiation"
            else (
                [0, 15]
                if meteo_var == "wind_speed"
                else ([0, 100] if meteo_var == "air_humidity" else [0, 40])
            )
        )
        axes[int(idx / nrows), idx % nrows].set_xlabel("")
        if meteo_var == "wind_speed":
            xticks = axes[int(idx / nrows), idx % nrows].get_xticks()
        axes[int(idx / nrows), idx % nrows].set_xticks(xticks)

        # axes[int(idx / nrows), idx % nrows].set_xticks([0, df.shape[0]])
        if idx < 2:
            axes[int(idx / nrows), idx % nrows].tick_params(length=0)
        axes[int(idx / nrows), idx % nrows].set_ylabel(
            meteo_vars[meteo_var], rotation=0, labelpad=13
        )
        axes[int(idx / nrows), idx % nrows].set_title(
            meteo_var.replace("_", " ").capitalize()
        )

    # _ = plt.xticks(rotation=0)
    plt.tight_layout()
    fig.set_size_inches(13, 6)
    fig.savefig(os.path.join(output_folder, "meteo.pdf"))
    fig.savefig(os.path.join(output_folder, "meteo.png"))


def water(
    tuning_folder=os.path.join("data", "errano_tuning"),
    evaluation_folder=os.path.join("data", "errano_evaluation_1gg"),
    output_folder=os.path.join("plots"),
):
    water_dict = {}
    water_vars = {
        "irrigation": "L",
        "precipitation": "mm",
    }

    for meteo_var in water_vars.keys():
        water_dict[meteo_var] = (
            pd.DataFrame({"timestamp": [1655251200], meteo_var: [0.0]})
            .append(
                pd.read_csv(os.path.join(tuning_folder, "water", f"{meteo_var}.csv")),
                ignore_index=True,
            )
            .append(
                pd.read_csv(
                    os.path.join(evaluation_folder, "water", f"{meteo_var}.csv")
                ),
                ignore_index=True,
            )
            .append(
                pd.DataFrame({"timestamp": [1661990400], meteo_var: [0.0]}),
                ignore_index=True,
            )
            .set_index("timestamp")
        )
    df = pd.concat(water_dict.values(), axis=1)
    df = df.reset_index()
    df["timestamp"] /= 3600 * 24
    df["timestamp"] = df["timestamp"].astype(int)
    df["timestamp"] *= 3600 * 24
    df["timestamp"] = pd.to_datetime(df["timestamp"], unit="s")
    df["timestamp"] = df["timestamp"].dt.strftime("%Y-%m-%d")
    df = df.rename(
        columns={
            "precipitation": "Precipitation in mm",
            "irrigation": "Irrigation in L",
        }
    )

    df = df.groupby("timestamp").sum()

    fig, ax = plt.subplots()
    df.plot(ax=ax, kind="bar")
    ax.set_xlabel("")
    ax.set_xticks([0, 15, 29, 46, 60, df.shape[0] - 1])
    ax.set_title("Precipitation and irrigation", fontdict={"fontsize": 18})
    ax2 = ax.twinx()
    ax.set_ylabel("L", rotation=0, labelpad=20, fontsize=15)
    ax2.set_ylabel("mm", rotation=0, labelpad=20, fontsize=15)
    ax2.set_ylim(0, 30)
    ax.set_ylim(0, 30)

    fig = ax.get_figure()
    fig.autofmt_xdate()

    plt.tight_layout()
    fig.set_size_inches(13, 6)
    fig.savefig(os.path.join(output_folder, "water.pdf"))
    fig.savefig(os.path.join(output_folder, "water.png"))


def ground_potential(
    tuning_folder=os.path.join("data", "errano_tuning"),
    evaluation_folder=os.path.join("data", "errano_evaluation_1gg"),
    output_folder=os.path.join("plots"),
):
    df = (
        pd.DataFrame({"timestamp": [1655251200]})
        .append(
            pd.read_csv(os.path.join(tuning_folder, "obs_data", "waterPotential.csv")),
            ignore_index=True,
        )
        .append(
            pd.read_csv(
                os.path.join(evaluation_folder, "obs_data", "waterPotential.csv")
            ),
            ignore_index=True,
        )
        .append(
            pd.DataFrame({"timestamp": [1661990400]}),
            ignore_index=True,
        )
    )

    df["timestamp"] = pd.to_datetime(df["timestamp"], unit="s")
    # df["timestamp"] = df["timestamp"].dt.strftime('%Y-%m-%d')
    df = df.set_index("timestamp")
    df = df.reindex(sorted(df.columns), axis=1)
    df = df.interpolate(method="linear", limit_direction="forward", axis=0)
    # print(df)

    # define subplot layout
    nrows, ncols = 3, 4
    fig, axes = plt.subplots(nrows=nrows, ncols=ncols)

    # add DataFrames to subplots

    for idx, meteo_var in reversed(list(enumerate(df.columns))):
        ax = axes[int(idx / ncols), idx % ncols]
        df[meteo_var].plot(ax=ax, sharex=True)
        ax.set_ylim([df.min().min(), -10])
        ax.set_yscale('symlog')
        ax.set_xlabel("")

        # ax.set_xticks([0,  df.shape[0] - 1])
        # if meteo_var == "z60_y0_x80":
        #     xticks = ax.get_xticks()
        # ax.set_xticks(xticks)

        # ax.set_xticks([0, df.shape[0]])
        if idx < 8:
            ax.tick_params(length=0)
        ax.set_ylabel("cbar")
        ax.set_title(
            meteo_var.replace("y0", "")
            .replace("_", " ")
            .replace("z", "depth = ")
            .replace("x", "distance = ")
        )

    # _ = plt.xticks(rotation=0)
    fig.set_size_inches(18, 6)
    plt.tight_layout()
    fig.savefig(os.path.join(output_folder, "ground_potential.pdf"))
    fig.savefig(os.path.join(output_folder, "ground_potential.png"))


def forecast_avg(
    obs_folder=os.path.join("data", "errano_evaluation_1gg"),
    forecast_folder=os.path.join("data", "errano_evaluation"),
    output_folder=os.path.join("plots"),
    with_forbidden_sensors=True
):
    support_dict = {
        "obs": os.path.join(obs_folder, "obs_data", "waterPotential.csv"),
    }
    for forecasting_day in ["1gg", "3gg", "7gg"]:
        support_dict[forecasting_day] = os.path.join(
            f"{forecast_folder}_{forecasting_day}", "output", "output.csv"
        )

    forecasting_dict = {}
    for data_type, input_path in support_dict.items():
        forecasting_dict[data_type] = (
            pd.DataFrame({"timestamp": [1655251200]})
            .append(
                pd.read_csv(input_path),
                ignore_index=True,
            )
            .set_index("timestamp")
            .drop(columns=[] if with_forbidden_sensors else forbidden_sensors)
        )
        forecasting_dict[data_type] *= -1
        forecasting_dict[data_type][forecasting_dict[data_type] < 20] = 20
        forecasting_dict[data_type] = forecasting_dict[data_type].apply(
            lambda x: np.log(x)
        )
        forecasting_dict[data_type] = forecasting_dict[data_type].reindex(
            sorted(forecasting_dict[data_type].columns), axis=1
        )
        forecasting_dict[data_type] = forecasting_dict[data_type].add_suffix(
            f"_{data_type}"
        )
    df = pd.concat(forecasting_dict.values(), axis=1)
    df = df.interpolate(method="linear", limit_direction="forward", axis=0)
    df = df.dropna(axis="index")
    df = df.reset_index()
    new_columns = []
    for data_type in support_dict.keys():
        if data_type != "obs":
            new_column = f"RMSE_{data_type}"
            new_columns += [new_column]
            # print(df[[c for c in df.columns if c.endswith("_obs")]].iloc[0])
            # print(df[[c for c in df.columns if c.endswith("_obs")]].iloc[:1])
            # print(df[[c for c in df.columns if c.endswith(f"_{data_type}")]].iloc[0])
            # print(df[[c for c in df.columns if c.endswith(f"_{data_type}")]].iloc[:1])
            df[new_column] = [
                mean_squared_error(
                    df[[c for c in df.columns if c.endswith("_obs")]].iloc[i:(i+1)],
                    df[[c for c in df.columns if c.endswith(f"_{data_type}")]].iloc[i:(i+1)],
                    squared=False,
                )
                for i in range(df.shape[0])
            ]
    df = df.append(pd.DataFrame({"timestamp": [1655251200]}), ignore_index=True)
    df = df.set_index("timestamp")
    df = df[new_columns]
    df.index = pd.to_datetime(df.index, unit="s")

    fig, ax = plt.subplots()
    df = df.rename(
        columns={
            column: column.replace("RMSE_", "forecasting horizon = ")
            for column in df.columns
        }
    )

    df.plot(ax=ax)
    ax.set_ylim([0, 2])
    ax.set_xlabel("")
    ax.set_ylabel("logRMSE")
    fig.set_size_inches(13, 6)
    plt.tight_layout()
    is_forbidden_sensors_string = "_with_forbidden_sensors" if with_forbidden_sensors else ""
    fig.savefig(os.path.join(output_folder, f"forecasting_avg{is_forbidden_sensors_string}.pdf"))
    fig.savefig(os.path.join(output_folder, f"forecasting_avg{is_forbidden_sensors_string}.png"))


def forecast_sensor(
    obs_folder=os.path.join("data", "errano_evaluation_1gg"),
    forecast_folder=os.path.join("data", "errano_evaluation"),
    output_folder=os.path.join("plots"),
    with_predicted=True
):
    support_dict = {
        "obs": os.path.join(obs_folder, "obs_data", "waterPotential.csv"),
    }
    for forecasting_day in ["1gg", "3gg", "7gg"]:
        support_dict[forecasting_day] = os.path.join(
            f"{forecast_folder}_{forecasting_day}", "output", "output.csv"
        )

    forecasting_dict = {}
    for data_type, input_path in support_dict.items():
        forecasting_dict[data_type] = (
            pd.DataFrame({"timestamp": [1655251200]})
            .append(
                pd.read_csv(input_path),
                ignore_index=True,
            )
            .set_index("timestamp")
            # .drop(columns=forbidden_sensors)
        )
        if not with_predicted:
            forecasting_dict[data_type] *= -1
            forecasting_dict[data_type][forecasting_dict[data_type] < 20] = 20
            forecasting_dict[data_type] = forecasting_dict[data_type].apply(
                lambda x: np.log(x)
            )
        else:
            forecasting_dict[data_type][forecasting_dict[data_type] > -20] = -20

        forecasting_dict[data_type] = forecasting_dict[data_type].reindex(
            sorted(forecasting_dict[data_type].columns), axis=1
        )
        forecasting_dict[data_type] = forecasting_dict[data_type].add_suffix(
            f"_{data_type}"
        )
    for forecasting_day in ["1gg", "3gg", "7gg"]:
        df = pd.concat([forecasting_dict["obs"], forecasting_dict[forecasting_day]], axis=1)
        df = df.interpolate(method="linear", limit_direction="forward", axis=0)
        df = df.dropna(axis="index")
        df = df.reset_index()
        sensor_columns = [c.replace("_obs", "") for c in df.columns if c.endswith("_obs")]
        if not with_predicted:
            new_columns = []
            for column in sensor_columns:
                new_column = f"RMSE_{column}"
                new_columns += [new_column]
                df[new_column] = [
                    mean_squared_error(
                        df[column + "_obs"].iloc[i:(i+1)],
                        df[column + f"_{forecasting_day}"].iloc[i:(i+1)],
                        squared=False,
                    )
                    for i in range(df.shape[0])
                ]
        df = df.append(pd.DataFrame({"timestamp": [1655251200]}), ignore_index=True)
        df = df.set_index("timestamp")
        if not with_predicted:
            df = df[new_columns]
        df.index = pd.to_datetime(df.index, unit="s")

        # define subplot layout
        nrows, ncols = 3, 4
        fig, axes = plt.subplots(nrows=nrows, ncols=ncols)

        for idx, meteo_var in reversed(list(enumerate(sensor_columns if with_predicted else df.columns))):
            ax = axes[int(idx / ncols), idx % ncols]
            columns = [meteo_var + suffix for suffix in ["_obs", f"_{forecasting_day}"]] if with_predicted else meteo_var
            df[columns].plot(ax=ax, sharex=True)
            ax.set_xlabel("")
            if not with_predicted:
                ax.set_ylim([0, 4])
                ax.set_ylabel("logRMSE")
            else:
                ax.set_ylim([df.min().min(), -10])
                ax.set_yscale('symlog')
                ax.set_ylabel("cbar")

            if idx < 8:
                ax.tick_params(length=0)
            ax.set_title(
                meteo_var.replace("RMSE_", "")
                .replace("y0", "")
                .replace("_", " ")
                .replace("z", "depth = ")
                .replace("x", "distance = "),
                color="red" if meteo_var.replace("RMSE_", "") in forbidden_sensors else "black"
            )

        fig.set_size_inches(18, 6)
        plt.tight_layout()
        suffix = "_with_predicted" if with_predicted else ""
        fig.savefig(os.path.join(output_folder, f"forecasting_sensor_{forecasting_day}{suffix}.pdf"))
        fig.savefig(os.path.join(output_folder, f"forecasting_sensor_{forecasting_day}{suffix}.png"))


def main():
    meteo()
    water()
    ground_potential()
    forecast_avg()
    forecast_avg(with_forbidden_sensors=False)
    forecast_sensor()
    forecast_sensor(with_predicted=False)


if __name__ == "__main__":
    main()
