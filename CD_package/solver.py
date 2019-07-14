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
# import pdb


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


def continuation_time(experiment_name,simulations_folder,V,dt):
    with open(simulations_folder+'{}/time_list.csv'.format(experiment_name), 'r') as csvfile:
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


############################### Create Logfile #################################






def SimpleNewtonSolve(gfu,a,tol=1e-13,maxits=25):
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

def errormess(mass):
    print("ERROR: Minimumvalue exploading or mass negative")
    print("Mass: {}".format(mass))


def paramlist(path,experiment_name,mass,alpha,epsilon,delta,dt,meshsize,order,T,geometry,ini_data_str):
    paramdict = {'experiment_name': experiment_name,'mass': mass, \
                 'alpha': alpha, 'epsilon': epsilon, 'delta': delta,\
                 'dt': dt,'meshsize': meshsize,'order': order, 'T': T, 'geometry': geometry,'ini_data': ini_data_str}
    with open(path+'{}/parameter.csv'.format(experiment_name),'w') as csvfile:
        fieldnames = ['experiment_name','mass','alpha','epsilon','delta','dt','meshsize','order','T','geometry','ini_data']
        paramwriter = csv.DictWriter(csvfile, fieldnames = fieldnames)
        paramwriter.writeheader()
        paramwriter.writerow(paramdict)

    return paramdict


def time_list_writer(experiment_name,dt,timestep,filename,linfmin_l2,linfmax_l2,mass):
    time_list_dict = {'dt': dt, 'timestep': timestep,'filename':filename,'linfymin': round(linfmin_l2,3),'linfmax':  round(linfmax_l2,3),'mass1': round(mass[0],3),'mass2':round(mass[1],3)}
    with open(filedir+'{}/time_list.csv'.format(experiment_name),'a') as csvfile:
        fieldnames = ['dt','timestep','filename','linfymin','linfmax','mass1','mass2']
        listwriter = csv.DictWriter(csvfile,fieldnames = fieldnames)
        listwriter.writerow(time_list_dict)


def run(experiment_name,uvold,gfucalc,gfuL2,a,mesh,dt,nsteps,paramlisto,startingtimestep,visualoutput_solver):
    initialmass = (Integrate(uvold.components[0],mesh),Integrate(uvold.components[1],mesh))
    with TaskManager():
        for i in range(startingtimestep,nsteps):
            atime = time()
            uvold.vec.data = gfucalc.vec
            # uvold.Set(gfucalc)
            print('inside')
            SimpleNewtonSolve(gfucalc,a)
            mass = (Integrate(gfucalc.components[0],mesh),Integrate(gfucalc.components[1],mesh))
            gfuL2.Set(gfucalc.components[0])
            linfmin_l2 = min(gfuL2.vec)
            linfmax_l2 = max(gfuL2.vec)
            modulo_constant = 1
            if dt < 1e-5:
                modulo_constant = 100
            elif dt < 1e-9:
                modulo_constant = 1000

            if i % modulo_constant == 0:
                gfucalc.Save(filedir+'{}/{}_time_{}'.format(experiment_name,experiment_name,str(i*dt)))
                with open(filedir+'{}/last_time.txt'.format(experiment_name),"w") as f:
                    f.write(str(i*dt))
                with open(filedir+'{}/last_time_before.txt'.format(experiment_name),"w") as f:
                    f.write(str((i-modulo_constant)*dt))
                filename_saved = filedir+'{}/{}_time_{}'.format(experiment_name,experiment_name,str(i*dt))
                print(experiment_name,dt,i,filename_saved)
                time_list_writer(experiment_name,dt,i,filename_saved,linfmin_l2,linfmax_l2,mass)


            if  linfmin_l2 < 0:
                print('>>>>>>>>>>>>>>>>> Negative Value >>>>>>>>>>>>>>>>> {}'.format(linfmin_l2))
                break

            print("total mass = ", mass, "MinvalueL2: ", linfmin_l2, "MaxvalueL2: ", linfmax_l2)
            if mass[0] > initialmass[0]*2  or linfmin_l2 <0:
                endtime = i*dt
                errormess(mass)
                with open(filedir+'{}/errormessage.txt'.format(experiment_name),"w") as f:
                    f.write('''Error at timestep {}:\nMass {}
                               \n min: {} \n max: {}'''.format(i,mass,linfmin_l2,linfmax_l2))

                with open(filedir+'{}/error_endtime.txt'.format(experiment_name),"w") as f:
                    f.write(str(int(i)))

                break
            if visualoutput_solver == True:
                Redraw(True)

            btime = time()
            print("step took: ",str(btime-atime))
            print("Time: {}, time step: {} of {}".format(str(dt*i),i,nsteps))
            endtime = nsteps*dt

        with open(filedir+'{}/parameter.txt'.format(experiment_name),"w") as f:
            for param in paramlisto:
                f.write(str(param) + ': '+str(paramlisto[param])+'\n')


    return gfucalc, endtime
