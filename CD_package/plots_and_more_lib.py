from netgen.geom2d import unit_square
from ngsolve import *
from ngsolve.internal import visoptions
from netgen.geom2d import SplineGeometry
import csv
import matplotlib.pyplot as plt
import os.path
import shutil
import pandas as pd
import numpy as np
from matplotlib.ticker import FormatStrFormatter,ScalarFormatter
import os
from pathlib import Path
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm


def reorder_timesteps(experiment_name,pathdirect):
    dirlist = sorted(Path(pathdirect).iterdir(), key=lambda f: f.stat().st_mtime)
#    pattern ='simulations/{}/{}_time_{}'.format(experiment_name,experiment_name,old_timestepping_last_time)
    timelist = []
    print('experiment name: '+experiment_name)
    cutlength = len(pathdirect+'{}_time_'.format(experiment_name))
    # print(cutlength)
    # print(pathdirect+'{}_time_'.format(experiment_name))
    for i in dirlist:
        tempstring = str(i)
        # pdb.set_trace()
        if tempstring[cutlength-2:] != '':
            timelist.append([float(tempstring[cutlength-2:]),str(i)])


    return sorted(timelist)


def opencsvfile(name,pathdirect):
    with open(pathdirect+'parameter.csv', 'r') as csvfile:
        csv_reader = csv.DictReader(csvfile)
        paramdict = list(csv_reader)

    return paramdict[0]
