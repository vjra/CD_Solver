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
    # moving to the narrow
    # midpointr = (0.1,0.5)
    # midpointr = (0.13,0.5)
    # rect
# Mesh parameters
# rect
    midpointr = (0.2,0.2)
    radius = 1.2
    midpointc = (0,0)
    meshsize = expridata['meshsize']
    geometry = 'unitcircle'
    order = 2
    T = expridata['T']
    dt = expridata['dt']
    length = 0
    # Bumpfunction parameters
    # beforemass1 = 13*pi
    mass1 = 14*pi
    mass2 = 0
    mass3 = 0
    xi = 1/500
    xi2 = 0

    # Equation parameters
    expri_number = expridata['exprinum']
    alpha = expridata['alpha']
    epsilon = expridata['epsilon']
    delta = expridata['delta']
    experiment_name = 'LIGHT_Expri_NOEV_refined_{}_alpha_{}_epsilon_{}_delta_{}'.format(expri_number,alpha,epsilon,delta)
    return  midpointr, radius, midpointc, meshsize, geometry, order, T, dt, xi, xi2, mass1, mass2,mass3, expri_number, alpha, epsilon, delta, experiment_name,length

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
for i in range(0,1,1):
    tempexpri = {'exprinum': int('20{}'.format(i)),'alpha': 1, 'delta': 0.0003, 'epsilon': 1, 'dt': 1e-5,'T': 0.7, 'meshsize': 0.008}
    experiment_list.append(tempexpri)
    print(tempexpri)


sleep(2)
############################### Create Logfile #################################
log_file_name = log_file_creator(simulations_folder,experiment_list)


############################### experiment loop #################################
for experiment in experiment_list:
    midpointr, radius, midpointc, meshsize, geometry, order, T, dt, xi, xi2, mass1, mass2,mass3, expri_number, alpha, epsilon, delta, experiment_name, length = fill_in_param(experiment)
    try:
        os.mkdir(filedir+experiment_name)
    except Exception as FileExistsError:
        pass
    counter_attempts = 0

    while counter_attempts < 1:
        try:
            # mesh,V,u,v,gfuL2 = meshgeneration_n_spaces_rec_trap(meshsize,order,length)
            mesh,V,u,v,gfuL2 = meshgeneration_n_spaces_rec_rect_edgy(meshsize,order,length)
            cont_switch = folder_checker(experiment_name,path)
            print('Contswitch variable: '+str(cont_switch))
            if cont_switch == True:
                old_timestepping_last_time, uvold = continuation_time(experiment_name,path,V,dt)
                startingtimestep = int(old_timestepping_last_time/dt)+1
            elif cont_switch == False:
                #old expri was with radial_two without baselevel
                uvold = initial_data_radial_two_baselevel(V,xi,xi2,mass1,mass2,midpointr,midpointc)
                uvold.Save(filedir+'{}/{}_time_{}'.format(experiment_name,experiment_name,0))
                startingtimestep = 1
            else:
                print('something went wrong in cont_switching')

            mass = Integrate(uvold,mesh)
            gfuL2.Set(uvold.components[0])
            linfmin_l2 = min(gfuL2.vec)
            linfmax_l2 = max(gfuL2.vec)
            if counter_attempts == 0:
                time_list_dict = {'dt': dt, 'timestep': 0,'filename':path+simulations_folder+'/'+experiment_name+'/'+experiment_name+'_time_0','linfymin': round(linfmin_l2,3),'linfmax':  round(linfmax_l2,3),'mass1': round(mass[0],3),'mass2':round(mass[1],3)}
                with open(filedir+'{}/time_list.csv'.format(experiment_name),'w') as csvfile:
                    fieldnames = ['dt','timestep','filename','linfymin','linfmax','mass1','mass2']
                    listwriter = csv.DictWriter(csvfile, fieldnames = fieldnames)
                    listwriter.writerow(time_list_dict)

            nsteps = int(floor(T/dt))
            print(nsteps)
            print("total mass = ", mass)
            a = weak_formulation(uvold,V,u,v,dt,delta,epsilon,alpha)
            if visualoutput_solver == True:
                #Redraw()
                gfucalc = GridFunction(V,name="u_n")
                gfucalc.Set(uvold)
                Draw(gfucalc)
                visoptions.scalfunction="u_n:1"
                visoptions.vecfunction = "None"
                visoptions.scaledeform1 = 0.001
                visoptions.deformation = 0
            else:
                gfucalc = GridFunction(V,name="u_n")
                gfucalc.Set(uvold)

            SetNumThreads(8)
            paramdicto = paramlist(experiment_name,xi,xi2,mass,mass2,alpha,epsilon,delta,dt,meshsize,order,T,geometry,midpointr,midpointc)
            gfucalc,endtime =run(experiment_name,uvold,gfucalc,gfuL2,a,mesh,dt,nsteps,paramdicto,startingtimestep,visualoutput_solver)
            print('Expri {} done'.format(expri_number))
            del a,mesh,V,u,v,uvold,gfucalc
            # dt = 10**(-(5+counter_attempts))
            counter_attempts += 1
            print('Attempts: {}, dt: {}, previous endtime: {}'.format(counter_attempts,dt,endtime))
            if endtime >= T:
                print('Done, endtime reached.')
                break

        except Exception as e:
            with open(log_file_name,'a') as f:
                f.write(str(e)+';\n')
        finally:
            pass
