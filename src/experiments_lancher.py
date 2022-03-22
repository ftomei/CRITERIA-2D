import subprocess
import time
import os
import argparse

import pandas as pd
import numpy as np

from hyperopt import fmin, tpe, hp, STATUS_OK, SparkTrials, Trials, pyll
from hyperopt.base import miscs_update_idxs_vals
from hyperopt.pyll.base import dfs, as_apply
from hyperopt.pyll.stochastic import implicit_stochastic_symbols
from sklearn.metrics import mean_squared_error, mean_absolute_error


iteration = 0

def objective(params):
    global iteration
    
    dataPath = os.path.join("data", "errano")
    outputPath = os.path.join(dataPath, "output")
    outputFileName = f"output_{str(iteration)}.csv"
    outputFilePath = os.path.join(outputPath, outputFileName)

    cmd = 'python src/main.py -it {} -vg_alpha {} -vg_n {} -theta_s {} -ks {} -root_depth_max {}  -root_x_def {} -root_z_def {} -kc_max {}'.format(
            iteration,
            params["VG_alpha"],
            params["VG_n"],
            params["thetaS"],
            params["Ks"],
            params["rootDepthMax"],
            params["rootXDeformation"],
            params["rootZDeformation"],
            params["kcMax"])
    subprocess.call(cmd, shell=True)

    simulated_data = pd.read_csv(outputFilePath)
    simulated_data = simulated_data.set_index('timestamp')
    simulated_data *= -1
    simulated_data[simulated_data < 20] = 20
    simulated_data = simulated_data.apply(lambda x: np.log(x))
    
    original_data = pd.read_csv(os.path.join(dataPath, 'ground_truth.csv'))
    original_data = original_data.set_index('timestamp')
    original_data = original_data.loc[simulated_data.index]
    original_data *= -1
    original_data[original_data < 20] = 20
    original_data = original_data.apply(lambda x: np.log(x))

    total_rmse = mean_squared_error(simulated_data, original_data, squared=False)
    
    iteration += 1
    return {'loss': total_rmse, 'status': STATUS_OK}


def main(args):
    print('Start')
    start_time = time.time()
    dataPath = os.path.join("data", "errano")
    outputPath = os.path.join(dataPath, "output")
    space = {
        'VG_alpha': hp.uniform('VG_alpha', 1.0, 3.0),
        'VG_n': hp.uniform('VG_n', 1.1, 1.4),
        'thetaS': hp.uniform('thetaS', 0.3, 0.5),
        'Ks': hp.loguniform('Ks', -16.1, -12.4),
        'rootDepthMax': hp.uniform('rootDepthMax', 0.6, 1.0),
        'rootXDeformation': hp.uniform('rootXDeformation', 0.0, 1.0),
        'rootZDeformation': hp.uniform('rootZDeformation', 0.0, 1.0),
        'kcMax': hp.uniform('kcMax', 1.5, 2.5),
        }
    trials = Trials()
    best = fmin(fn=objective,
                space=space,
                algo=tpe.suggest,
                max_evals=args.num_iterations,
                trials=trials,
                show_progressbar=True,
                rstate=np.random.RandomState(42))

    best_parameters = pd.DataFrame.from_records([best])
    best_parameters.to_csv(os.path.join(outputPath, 'best_parameters.csv'), index=False)

    all_trials = pd.DataFrame(trials.trials)
    all_trials.to_json(os.path.join(outputPath, 'all_trials.json'), indent=True)

    best_trial = pd.DataFrame(trials.best_trial)
    best_trial.to_json(os.path.join(outputPath,'best_trial.json'), indent=True)

    end_time = time.time()
    print("Time of execution:", end_time-start_time)

def parse_args():

    parser = argparse.ArgumentParser(description="CRITERIA-3D - Soil simulator calibration")

    # I need to force argparse to do what was descriped in spec

    parser.add_argument("-nits", "--num_iterations", nargs="?", type=int, required=True,
                        help="number of iterations")

    args = parser.parse_args()

    return args

if __name__ == "__main__":
    args = parse_args()
    main(args)