from netgen.geom2d import unit_square
from netgen.geom2d import SplineGeometry
from ngsolve import *
from math import pi
from ngsolve.internal import visoptions
from time import time
import os
import csv
import numpy as np
from netgen.meshing import PointId
from time import sleep
from datetime import datetime
from pathlib import Path
import fnmatch
import logging

import pdb

############################################################################
########################## New features, not implemented ###################
############################################################################

############################################################################
########################## Logging #########################################
############################################################################

# create and configure logger
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)

# create a file handler
handler = logging.FileHandler('./logs/solver.log', mode='w')
handler.setLevel(logging.DEBUG)

# create a logging format
formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
handler.setFormatter(formatter)

# add the file handler to the logger
logger.addHandler(handler)

############################################################################
########################## Functions to log the experiment #################
############################################################################

def folder_checker(experiment_name,path):
    """Checks if a folder for the experiment (experiment_name)
     exists. If this is true, it sets the continuation parameter
     cont_switch to 'True'. If this false, creates a
     folder for the experiment in the directory 'simulations_data'.

    Parameters:
    experiment_name (str): url that contains a pdf.
    path (str): path to the folder for the simulation data.

    Returns:
    cont_switch (boolean): is false if no experiment folder is present,
                           and True if its a new experiment attempt.
    """

    cont_switch=False
    uvoldname = ''
    try:
        os.mkdir(path+'{}'.format(experiment_name,experiment_name))

    except Exception as FileExistsError:
        try:
            if len(os.listdir(path+'/{}'.format(experiment_name))) == 0:
                cont_switch = False
            else:
                cont_switch = True
        except OSError as e:
            print ("Error: %s - %s." % (e.filename, e.strerror))
    return cont_switch


def continuation_time(experiment_name,path,V,dt,uvold):
    """Extracts the stored data for the continuation of an experiment.

    Parameters:
    experiment_name (str): name of the experiment to load.
    path (str): path to the folder for the simulation data.
    V (ngsolve.comp.FESpace object): FEspace for this experiment.
    dt (float): time step size for this experiment.
    uvold (ngsolve.comp.GridFunction object): gridfunction of last time step stored in folder.

    Returns:
    old_timestepping_last_time (float): last time step stored in folder.
    uvold (ngsolve.comp.GridFunction object): gridfunction of last time step stored in folder.
    """
    with open(path+'{}/time_list.csv'.format(experiment_name), 'r') as csvfile:
        reader = csv.reader(csvfile)
        lines = list(reader)

        if 0 < len(lines) <= 2:
            row = lines[-1]
            old_timestepping_last_time = 0
        elif len(lines)>2:
            row = lines[-2]
            old_timestepping_last_time = float(row[1])*float(row[0])
        else:
            old_timestepping_last_time = 0
            row = []

    print('We continue according the the data:')
    print(row)
    print('We continue example at time step: {}.\n'.format(old_timestepping_last_time))
    if len(lines) > 1:
        uvold.Load(str(row[2]))

    return old_timestepping_last_time, uvold

def log_file_creator(path,experiment_list):
    """Creates a log file, with date and time when the experiment was started.
    Contains the experiment name and error if occured.

    Parameters:
    experiment_name (str): url that contains a pdf.
    path (str): path to the folder for the simulation data.

    Returns:
    logfilename (str): name of the logfile.
    """

    today = datetime.now()
    today_string = '{}_{}_{}_{}'.format(today.year,today.month,today.day,today.hour)
    logfilename = path+'{}_errorlog.csv'.format(today_string)
    try:
        with open(logfilename,"w") as f:
            f.write("##############LOG############\n")
            experiment_string =''
            for key,value in experiment_list.items():
                experiment_string += str(key) +': '+ str(value) +', '

            f.write(experiment_string+';\n')
    except Exception as FileExistsError:
        print(FileExistsError)

    return logfilename

def paramlist(path,experiment_name,mass,alpha,epsilon,delta,dt,meshsize,order,T,geometry,ini_data_str = ''):
    """Creates a parameter dictionary that contains all parameter set.
    In addition, a csv file (parameter.csv) is created in the experiment folder,
    containing the same information.

    Parameters:
    experiment_name (str): url that contains a pdf.
    path (str): path to the folder for the simulation data.
    mass (float): L^1 norm of the initial data.
    dt (float): time step size.
    meshsize (float): mesh size.
    order (float): order of the functions in the finite element space.
    T (float): end time of the simulation.
    geometry (str): geometry specification as string.
    ini_data_str (str): initial data as string, if given.

    alpha, epsilon, delta (str): parameters for the system in equation, i.e.
        rho_t         = div(\nabla \rho - rho*\nabla c)
        epsilon*c_t   = div(\nabla c + delta*\nabla rho ) + c + rho^alpha


    Generates:
    parameter.csv: file containing all parameters of the experiment.
    Returns:
    paramdict (dict): dictionary containing all parameters of the experiment.
    """
    paramdict = {'experiment_name': experiment_name,'mass': mass, \
                 'alpha': alpha, 'epsilon': epsilon, 'delta': delta,\
                 'dt': dt,'meshsize': meshsize,'order': order, 'T': T, 'geometry': geometry,'ini_data': ini_data_str}
    with open(path+'{}/parameter.csv'.format(experiment_name),'w') as csvfile:
        fieldnames = ['experiment_name','mass','alpha','epsilon','delta','dt','meshsize','order','T','geometry','ini_data']
        paramwriter = csv.DictWriter(csvfile, fieldnames = fieldnames)
        paramwriter.writeheader()
        paramwriter.writerow(paramdict)

    return paramdict


def time_list_writer(path,experiment_full_name,dt,timestep,filename,linfmin_l2,linfmax_l2,mass):
    """Creates a timelist to restart the experiment at the last time step recorded.

    Parameters:
    experiment_full_name (str): url that contains a pdf.
    path (str): path to the folder for the simulation data.
    dt (float): time step size.
    timestep (int): current time step.
    filename (str): file name of the time step computed.
    linfmin_l2 (float): approximation of the minimum of rho.
    linfmax_l2 (float): approximation of the maximum of rho.
    mass (float): L^1 norm of the initial data.

    Generates:
    time_list.csv: list of all time step taken, including some information like mass.

    """
    time_list_dict = {'dt': dt, 'timestep': timestep,'filename':filename,'linfymin': round(linfmin_l2,3),'linfmax':  round(linfmax_l2,3),'mass1': round(mass[0],3),'mass2':round(mass[1],3)}
    with open(path+'{}/time_list.csv'.format(experiment_full_name),'a') as csvfile:
        print(path+'{}/time_list.csv'.format(experiment_full_name))
        fieldnames = ['dt','timestep','filename','linfymin','linfmax','mass1','mass2']
        listwriter = csv.DictWriter(csvfile,fieldnames = fieldnames)
        listwriter.writerow(time_list_dict)

############################################################################
########################## Solver functions ################################
############################################################################

def SimpleNewtonSolve(gfu,a,tol=1e-13,maxits=25):
    """Newton solver for nonlinear systems, uses the ngsolve routine
    'AssembleLinearization' to compute a linearization.

    Parameters:
    gfu (ngsolve.comp.GridFunction): solution from last time step.
    a (ngsolve.comp.BilinearForm): bilinear form of the divergence term.
    tol (float): error tolerance
    maxits (int): maximum of iterations before the solver stops in any case.

    Updates:
    gfu.vec.data (ngsolve.comp.GridFunction)

    """
    res = gfu.vec.CreateVector()
    du = gfu.vec.CreateVector()
    fes = gfu.space
    for it in range(maxits):
        print ("Iteration {:3}  ".format(it),end="")
        a.Apply(gfu.vec, res)
        a.AssembleLinearization(gfu.vec)
        du.data = a.mat.Inverse(fes.FreeDofs()) * res
        gfu.vec.data -= du

        #stopping criteria
        stopcritval = sqrt(abs(InnerProduct(du,res)))
        print ("<A u",it,", A u",it,">_{-1}^0.5 = ", stopcritval)
        if stopcritval < tol:
            break





def run(path,experiment_full_name,uvold,gfucalc,gfuL2,a, \
        mesh,dt,nsteps,paramlisto,startingtimestep,visualoutput_solver = False,log_file_name = ''):
    """Actual solver routine, that loops over the time steps and uses
    Newton in the function 'SimpleNewtonSolve' to solve the nonlinear system.

    Parameters:
    experiment_full_name (str): url that contains a pdf.
    path (str): path to the folder for the simulation data.
    uvold (ngsolve.comp.GridFunction): initial data as grid function.
    gfucalc (ngsolve.comp.GridFunction): grid function on FEspace that stores latest time step of the solution.
    gfuL2 (ngsolve.comp.GridFunction): grid function on L2 space that stores latest time step of the solution.
    a (ngsolve.comp.BilinearForm): bilinear form of the divergence term.
    mesh (ngsolve.comp.Mesh): mesh.
    dt (float): time step size.
    nsteps (int): number of time step to take.
    paramlisto (dict): dictionary of parameters, see function 'paramlist'.
    startingtimestep (int): initial time step.
    visualoutput_solver (boolean): switch to activate plotting in netgen.
    log_file_name (str): filename of log file.

    Return:
    gfucalc (ngsolve.comp.GridFunction): grid function of solution at last time step on FEspace.
    endtime (float): time when routine stopped.
    """


    initialmass = (Integrate(uvold.components[0],mesh),Integrate(uvold.components[1],mesh))
    with TaskManager():
        for i in range(startingtimestep,nsteps):
            atime = time()
            endtime = str(dt*i)
            uvold.vec.data = gfucalc.vec
            SimpleNewtonSolve(gfucalc,a)
            mass = (Integrate(gfucalc.components[0],mesh),Integrate(gfucalc.components[1],mesh))
            gfuL2.Set(gfucalc.components[0])
            linfmin_l2 = min(gfuL2.vec)
            linfmax_l2 = max(gfuL2.vec)
            # todo: implement modulo constant, in order to only save every time step modulo this constant.
            modulo_constant = 0
            gfucalc.Save(path+'{}/{}_time_{}'.format(experiment_full_name,experiment_full_name,str(i*dt)))
            with open(path+'{}/last_time.txt'.format(experiment_full_name),"w") as f:
                f.write(str(i*dt))
            with open(path+'{}/last_time_before.txt'.format(experiment_full_name),"w") as f:
                f.write(str((i-modulo_constant)*dt))
            filename_saved = path+'{}/{}_time_{}'.format(experiment_full_name,experiment_full_name,str(i*dt))
            print('dt: ',dt,'filename:',filename_saved)
            time_list_writer(path,experiment_full_name,dt,i,filename_saved,linfmin_l2,linfmax_l2,mass)
            print("total mass = ", mass, "MinvalueL2: ", linfmin_l2, "MaxvalueL2: ", linfmax_l2)
            if linfmin_l2 < 0:
                print('>>>>>>>>>>>>>>>>> Negative Value >>>>>>>>>>>>>>>>> {}'.format(linfmin_l2))
                with open(log_file_name,'a') as f:
                    f.write(""">>>>>>>>>>>>>>>>> Reached negative value in rho
                            (gfucalc[0]) >>>>>>>>>>>>>>>>> {}""".format(linfmin_l2)+';\n')

                logger.info('''Error at timestep {}:\nMass {}
                           \n min: {} \n max: {}'''.format(i,mass,linfmin_l2,linfmax_l2))


                endtime = i*dt
                with open(path+'{}/errormessage.txt'.format(experiment_full_name),"w") as f:
                    f.write('''Error at timestep {}:\nMass {}
                               \n min: {} \n max: {}'''.format(i,mass,linfmin_l2,linfmax_l2))

                with open(path+'{}/error_endtime.txt'.format(experiment_full_name),"w") as f:
                    f.write(str(int(i)))
                raise Exception('Negative value reached.')

            if visualoutput_solver == True:
                Redraw(True)

            btime = time()
            print("step took: ",str(btime-atime))
            print("Time: {}, time step: {} of {}".format(str(dt*i),i,nsteps))


        with open(path+'{}/parameter.txt'.format(experiment_full_name),"w") as f:
            for param in paramlisto:
                f.write(str(param) + ': '+str(paramlisto[param])+'\n')


    return gfucalc, endtime
