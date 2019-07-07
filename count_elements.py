from netgen.geom2d import unit_square
from ngsolve import *
from math import pi
from ngsolve.internal import visoptions
from time import time
from time import sleep
from datetime import datetime
from ks_solver4_new_many_v2 import *
import os
import csv
from pathlib import Path
import fnmatch

############################### auxilary functions #################################

path = './'
visualoutput_solver = True
mesh_switch = 'edgy'

if visualoutput_solver == True:
    import netgen.gui

def folder_checker(experiment_name,path):
    cont_switch = False
    uvoldname = ''
    try:
        os.mkdir(path+'simulations/{}'.format(experiment_name,experiment_name))

    except Exception as FileExistsError:
        try:
            if len(os.listdir(path+'simulations/{}'.format(experiment_name))) == 0:
                cont_switch = False
            else:
                cont_switch = True
        except OSError as e:
            print ("Error: %s - %s." % (e.filename, e.strerror))
    return cont_switch

def continuation_time(experiment_name,path,V,dt):
    with open(path+'simulations/{}/time_list.csv'.format(experiment_name), 'r') as csvfile:
        reader = csv.reader(csvfile)
        lines = list(reader)
        if len(lines) == 2:
            row = lines[-1]
            old_timestepping_last_time = 0
        else:
            row = lines[-2]
            old_timestepping_last_time = float(row[1])*float(row[0])

    print('We continue according the the data:')
    print(row)
    uvold = GridFunction(V)
    print('We continue example {}. \n At timestep:'.format(old_timestepping_last_time))
    print(str(old_timestepping_last_time))
    uvold.Load(str(row[2]))


    return old_timestepping_last_time, uvold



def fill_in_param(expridata):
    # Mesh parameters
    midpointr = (0.5,0.5)
    radius = 1
    midpointc = (0,0)
    meshsize = expridata['meshsize']
    geometry = 'unitcircle'
    order = 2
    T = expridata['T']
    dt = expridata['dt']

    # Bumpfunction parameters
    xi = 1e-2
    xi2 = 1
    mass = 9*pi
    mass2 = 0

    # Equation parameters
    expri_number = expridata['exprinum']
    alpha = expridata['alpha']
    epsilon = expridata['epsilon']
    delta = expridata['delta']
    experiment_name = 'Hitt_Expri_{}_alpha_{}_epsilon_{}_delta_{}'.format(expri_number,alpha,epsilon,delta)
    return midpointr, radius, midpointc, meshsize, geometry, order, T, dt, xi, xi2, mass, mass2, expri_number, alpha, epsilon, delta, experiment_name

def log_file_creator(foldername,experiment_list):
    today = datetime.now()
    today_string = '{}_{}_{}_{}'.format(today.year,today.month,today.day,today.hour)
    logfilename = foldername+'/{}_errorlog.csv'.format(today_string)
    try:
        with open(logfilename,"w") as f:
            f.write("##############LOG############\n")
            for expridata in experiment_list:
                f.write(str(expridata)+';\n')
    except Exception as FileExistsError:
        pass

    return logfilename

############################### global variables #################################
# subfolder where the simulations are stored
simulations_folder = 'simulations'

# !!!! visualoutput_solver is a global variable in the module ks_solver!!!!

############################### experiments parameter #################################
experiment_list = []
for i in range(2,5,1):
    tempexpri = {'exprinum': int('25{}'.format(i)),'alpha': 2.5, 'delta': 10**(-(i)) , 'epsilon': 1, 'dt': 1e-2,'T': 5, 'meshsize': 0.035}
    experiment_list.append(tempexpri)
    print(tempexpri)


sleep(2)
############################### Create Logfile #################################
log_file_name = log_file_creator(simulations_folder,experiment_list)


############################### experiment loop #################################
for experiment in experiment_list:
    midpointr, radius, midpointc, meshsize, geometry, order, T, dt, xi, xi2, mass, mass2, expri_number, alpha, epsilon, delta, experiment_name = fill_in_param(experiment)
    try:
        os.mkdir(filedir+experiment_name)
    except Exception as FileExistsError:
        pass

    if mesh_switch == 'edgy':
        mesh,V,u,v,gfuL2 = meshgeneration_n_spaces_rec_edgy(meshsize,order,midpointc,radius)
    elif mesh_switch == 'center':
        mesh,V,u,v,gfuL2 = meshgeneration_n_spaces_rec_center(meshsize,order,midpointc,radius)
    else:
        mesh,V,u,v,gfuL2 = meshgeneration_n_spaces_rec(meshsize,order,midpointc,radius)
