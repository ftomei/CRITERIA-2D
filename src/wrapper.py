import subprocess
import time
import os
import json
import psutil


import pandas as pd
import numpy as np

from functools import partial
from flaml import tune
from sklearn.metrics import mean_squared_error

from tuning.space_loading import get_space
from tuning.utils import create_directory, parse_args_tuning
from tuning.json_to_csv import json_to_csv

iteration = 1


def kill(proc_pid):
    process = psutil.Process(proc_pid)
    for proc in process.children(recursive=True):
        proc.kill()
    process.kill()


def objective(dataPath, params):
    global iteration
    result = {"rmse": float("inf"), "status": "fail"}

    try:
        iterations_str = f"_{iteration}" if params else ""
        outputPath = os.path.join(dataPath, "output")
        inputFilePath = os.path.join(outputPath, f"input{iterations_str}.json")
        outputFilePath = os.path.join(outputPath, f"output{iterations_str}.csv")
        stdoutFilePath = os.path.join(outputPath, f"stdout{iterations_str}.txt")
        stderrFilePath = os.path.join(outputPath, f"stderr{iterations_str}.txt")
        observedDataPath = os.path.join(dataPath, "obs_data", "waterPotential.csv")

        if params:
            json.dump(params, open(inputFilePath, "w"))

        open(stdoutFilePath, "w")
        open(stderrFilePath, "w")
        with open(stdoutFilePath, "a") as log_out:
            with open(stderrFilePath, "a") as log_err:
                try:
                    process = subprocess.Popen(
                        f"python src/main.py -p {dataPath} -it {iteration if params else -1}",
                        shell=True,
                        stdout=log_out,
                        stderr=log_err,
                    )
                    process.wait(timeout=5400)
                except Exception as e:
                    # print(e)
                    kill(process.pid)

        simulated_data = pd.read_csv(outputFilePath)
        simulated_data = simulated_data.set_index("timestamp")
        simulated_data = simulated_data.drop(columns=["z20_y0_x0", "z20_y0_x25", "z40_y0_x0", "z40_y0_x50"])
        simulated_data *= -1
        simulated_data[simulated_data < 20] = 20
        simulated_data = simulated_data.apply(lambda x: np.log(x))

        original_data = pd.read_csv(observedDataPath)
        original_data = original_data.set_index("timestamp")
        original_data = original_data.drop(columns=["z20_y0_x0", "z20_y0_x25", "z40_y0_x0", "z40_y0_x50"])
        original_data = original_data.loc[simulated_data.index]
        original_data *= -1
        original_data[original_data < 20] = 20
        original_data = original_data.apply(lambda x: np.log(x))

        result["rmse"] = mean_squared_error(
            simulated_data, original_data, squared=False
        )
        result["status"] = "success"
    except Exception as e:
        print(
            f"""MyException: {e}"""
            #   {traceback.print_exc()}"""
        )

    iteration += 1
    return result


def main(args):
    print("Start")
    np.random.seed(args.seed)
    settingsFolder = os.path.join(args.path, "settings")
    outputFolder = create_directory(args.path, "output")
    outputFolder = create_directory(outputFolder, "summary")
    outputFile = os.path.join(outputFolder, "summary")
    space = get_space(os.path.join(settingsFolder, "tuning_space.json"))

    start_time = time.time()
    if args.num_iterations <= 0:
        result = objective(args.path, params=None)
    else:
        analysis = tune.run(
            evaluation_function=partial(objective, args.path),
            config=space,
            metric="rmse",
            mode="min",
            num_samples=args.num_iterations,
            verbose=1,
            max_failure=args.num_iterations,
        )
    end_time = time.time()

    outcome = [result] if args.num_iterations <= 0 else analysis.results.values()

    # Specify which information are needed for the output
    filtered_keys = ["rmse", "status", "config", "time_total_s"]
    # Prepare the output file
    automl_output = {
        "optimization_time": end_time - start_time,
        # For each visited config, filter the information
        "results": [
            {
                key: value if value != float("-inf") else str(value)
                for key, value in values.items()
                if key in filtered_keys
            }
            for values in outcome
        ],
    }

    if args.num_iterations > 0:
        automl_output["best_config"] = {
            key: value
            for key, value in analysis.best_trial.last_result.items()
            if key in filtered_keys
        }

    # Export the result
    with open(f"{outputFile}.json", "w") as outfile:
        json.dump(automl_output, outfile)

    # Convert the result in csv
    json_to_csv(automl_output=automl_output.copy(), outputFile=f"{outputFile}.csv")


args = parse_args_tuning()
main(args)
