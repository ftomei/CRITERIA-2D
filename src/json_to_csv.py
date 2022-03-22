import os
import json
import argparse

import pandas as pd
import numpy as np

def main(args):
    dataPath = os.path.join("data", "errano")
    outputPath = os.path.join(dataPath, "output")
    allTrials = os.path.join(outputPath, "all_trials.json")
    df = pd.DataFrame()
    with open(allTrials) as json_file:
        data = json.load(json_file)
        for i in range(args.num_iterations):
            df = df.append({
                'loss': data['result'][str(i)]['loss'],
                'Ks': data['misc'][str(i)]['vals']['Ks'][0],
                'VG_alpha': data['misc'][str(i)]['vals']['VG_alpha'][0],
                'VG_n': data['misc'][str(i)]['vals']['VG_n'][0],
                'kcMax': data['misc'][str(i)]['vals']['kcMax'][0],
                'rootDepthMax': data['misc'][str(i)]['vals']['rootDepthMax'][0],
                'rootXDeformation': data['misc'][str(i)]['vals']['rootXDeformation'][0],
                'rootZDeformation': data['misc'][str(i)]['vals']['rootZDeformation'][0],
                'thetaS': data['misc'][str(i)]['vals']['thetaS'][0],
                }, ignore_index=True)
    df.insert(0, 'loss', df.pop('loss'))
    df.to_csv(os.path.join(outputPath, 'result.csv'))

def parse_args():

    parser = argparse.ArgumentParser(description="CRITERIA-3D - Results collection")

    # I need to force argparse to do what was descriped in spec

    parser.add_argument("-nits", "--num_iterations", nargs="?", type=int, required=True,
                        help="number of iterations")

    args = parser.parse_args()

    return args

if __name__ == "__main__":
    args = parse_args()
    main(args)