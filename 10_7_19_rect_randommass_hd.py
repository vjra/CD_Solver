from netgen.geom2d import unit_square
from ngsolve import *
from math import pi
from ngsolve.internal import visoptions
from time import time
from time import sleep
from datetime import datetime
from ks_solver5_v1 import *
import os
import csv
from pathlib import Path
import fnmatch


############################### auxilary functions #################################







def fill_in_param(expridata):
    midpointr = (0,0)
    radius = 1
    midpointc = (0,0)
    meshsize = expridata['meshsize']
    geometry = 'unitcircle'
    order = 2
    T = expridata['T']
    dt = expridata['dt']
    length = 0
    # Bumpfunction parameters
    mass1 = 0.05
    n_of_drops = 2000
    xi = 1/800
    # Equation parameters
    expri_number = expridata['exprinum']
    alpha = expridata['alpha']
    epsilon = expridata['epsilon']
    delta = expridata['delta']
    experiment_name = 'RECTMASS_Expri_HD_{}_alpha_{}_epsilon_{}_delta_{}'.format(expri_number,alpha,epsilon,delta)
    return  midpointr, radius, midpointc, meshsize, geometry, order, T, dt, xi, mass1,n_of_drops, expri_number, alpha, epsilon, delta, experiment_name,length











############################### experiment loop #################################
midpointr, radius, midpointc, meshsize, geometry, order, T, dt, xi, mass1,n_of_drops, expri_number, alpha, epsilon, delta, experiment_name, length = fill_in_param(experiment)
try:
    os.mkdir(filedir+experiment_name)
except Exception as FileExistsError:
    pass
counter_attempts = 0

mesh,V,u,v,gfuL2 = meshgeneration_n_spaces_rect_hd_c_refine(meshsize,order)

# mesh,V,u,v,gfuL2 = meshgeneration_n_spaces_rect_hd_c(meshsize,order)
cont_switch = folder_checker(experiment_name,path)
print('Contswitch variable: '+str(cont_switch))
if cont_switch == True:
    old_timestepping_last_time, uvold = continuation_time(experiment_name,path,V,dt)
    startingtimestep = int(old_timestepping_last_time/dt)+1

elif cont_switch == False:
    uvold =  initial_data_random_many(mass1,V,xi,n_of_drops)
    uvold.Save(filedir+'{}/{}_time_{}'.format(experiment_name,experiment_name,0))
    startingtimestep = 1
else:
    print('something went wrong in cont_switching')

mass = Integrate(uvold.components[0],mesh)
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
    gfucalc.vec.data = uvold.vec
    Draw(gfucalc.components[0])
    # visoptions.scalfunction="u_n:1"
    visoptions.vecfunction = "None"
    visoptions.scaledeform1 = 0.001
    visoptions.deformation = 0
else:
    gfucalc = GridFunction(V,name="u_n")
    #gfucalc.components[0].Set(uvold.components[0])
    gfucalc.vec.data = uvold.vec

SetNumThreads(8)
paramdicto = paramlist(experiment_name,xi,mass,n_of_drops,alpha,epsilon,delta,dt,meshsize,order,T,geometry,midpointr,midpointc)

gfucalc,endtime =run(experiment_name,uvold,gfucalc,gfuL2,a,mesh,dt,nsteps,paramdicto,startingtimestep,visualoutput_solver)
print('Expri {} done'.format(expri_number))
del a,mesh,V,u,v,uvold,gfucalc
