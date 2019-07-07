from netgen.geom2d import unit_square
from ngsolve import *
from math import pi
from ngsolve.internal import visoptions
from time import time
from time import sleep
from datetime import datetime
from ks_solver4_new_EV import *
import os
import csv
from pathlib import Path
import fnmatch
import sys

############################### auxilary functions #################################

path = './'
visualoutput_solver = True
mesh_switch = 'rect'

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
        row = lines[-1]
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
    # rect
    midpointr = (0.2,0.2)
    radius = 1.2
    midpointc = (0,0)
    meshsize = expridata['meshsize']
    geometry = 'unitcircle'
    order = 3
    T = expridata['T']
    dt = expridata['dt']
    length = 0
    # Bumpfunction parameters
    mass1 = 13*pi
    mass2 = 0
    mass3 = 0
    xi = 1/500
    xi2 = 0

    # Equation parameters
    expri_number = expridata['exprinum']
    alpha = expridata['alpha']
    epsilon = expridata['epsilon']
    delta = expridata['delta']
    experiment_name = 'LIGHT_Expri_{}_alpha_{}_epsilon_{}_delta_{}'.format(expri_number,alpha,epsilon,delta)
    return midpointr, radius, midpointc, meshsize, geometry, order, T, dt, xi, xi2, mass1, mass2,mass3, expri_number, alpha, epsilon, delta, experiment_name,length

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
for i in range(5,6,1):
    tempexpri = {'exprinum': int('4{}'.format(i)),'alpha': 1, 'delta': 0.0001 , 'epsilon': 1, 'dt': 1e-4,'T': 1, 'meshsize': 0.005}
    experiment_list.append(tempexpri)
    print(tempexpri)


sleep(2)
############################### Create Logfile #################################
log_file_name = log_file_creator(simulations_folder,experiment_list)

print(experiment_list)
#
############################### experiment loop #################################
for experiment in experiment_list:
    midpointr, radius, midpointc, meshsize, geometry, order, T, dt, xi, xi2, mass1, mass2,mass3, expri_number, alpha, epsilon, delta, experiment_name,length = fill_in_param(experiment)
    try:
        os.mkdir(filedir+experiment_name)
    except Exception as FileExistsError:
        pass

    try:
        if mesh_switch == 'edgy':
            mesh,V,u,v,gfuL2 = meshgeneration_n_spaces_rec_edgy(meshsize,order,midpointc,radius)
        elif mesh_switch == 'center':
            mesh,V,u,v,gfuL2 = meshgeneration_n_spaces_rec_center(meshsize,order,midpointc,radius)
        elif mesh_switch == 'obst':
            mesh,V,u,v,gfuL2 = meshgeneration_n_spaces_rec_obst(meshsize,order,midpointc,radius)
        elif mesh_switch == 'rect':
            mesh,V,u,v,gfuL2 = meshgeneration_n_spaces_rec_rect(meshsize,order,length)
        else:
            mesh,V,u,v,gfuL2 = meshgeneration_n_spaces_rec(meshsize,order,midpointc,radius)

        cont_switch = folder_checker(experiment_name,path)
        print('Contswitch variable: '+str(cont_switch))
        if cont_switch == True:
            old_timestepping_last_time, rhovold = continuation_time(experiment_name,path,V,dt)
            startingtimestep = int(old_timestepping_last_time/dt)+1
        elif cont_switch == False:
            print('yes')
            rhovold =  initial_data_radial_two(V,xi,xi2,mass1,mass2,midpointr,midpointc)
            rhovold.Save(filedir+'{}/{}_time_{}'.format(experiment_name,experiment_name,0))
            startingtimestep = 1
        else:
            print('something went wrong in cont_switching')

        gfucalcw = GridFunction(V)
        gfucalcw.Set(CoefficientFunction((log(rhovold[0])*delta,rhovold[1])))
        gfucalc = GridFunction(V,name="rhov")
        gfucalc.Set(rhovold)
        mass = Integrate(rhovold,mesh)
        nsteps = int(ceil(T/dt))+1
        print(nsteps)
        print("total mass = ", mass)
        b = weak_formulation_EV(rhovold,V,u,v,dt,delta,epsilon,alpha)
        if visualoutput_solver == True:
            Draw(gfucalc)
            visoptions.scalfunction="rhov:1"
            visoptions.vecfunction = "None"
            visoptions.scaledeform1 = 0.1
            visoptions.deformation = 0
        else:
            gfucalc = GridFunction(V,name="rhov")
            gfucalc.Set(rhovold)

        SetNumThreads(8)
        paramdicto = paramlist(experiment_name,xi,xi2,mass,mass2,alpha,epsilon,delta,dt,meshsize,order,T,geometry,midpointr,midpointc)
        atime = time()
        gfucalc,endtime =run_EV(experiment_name,rhovold,gfucalc,gfucalcw,gfuL2,b,mesh,dt,nsteps,delta,paramdicto,startingtimestep,visualoutput_solver)
        btime = time()
        print('Expri {} done, total time: {}'.format(expri_number,totaltime))
        with open(path+'simulations/{}/parameter.txt'.format(experiment_name),'a') as f:
            f.write('total time: {}'.format(totaltime))
        del b,mesh,V,u,v,uvold,gfucalc,gfucalcw

    except Exception as e:
        with open(log_file_name,'a') as f:
            f.write('Error on line {}'.format(sys.exc_info()[-1].tb_lineno)+'; '+str(type(e).__name__)+'; ' +str(e)+';\n')
    finally:
        pass
