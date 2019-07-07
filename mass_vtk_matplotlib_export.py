from netgen.geom2d import unit_square
from ngsolve import *
from ngsolve.internal import visoptions
from netgen.geom2d import SplineGeometry
# from ks_solver4_new_EV import *
# from ks_solver4_new_many_v2 import *
from ks_solver5_v1 import *
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
    print('expriname: '+experiment_name)
    cutlength = len(pathdirect+'{}_time_'.format(experiment_name))
    print(cutlength)
    print(pathdirect+'{}_time_'.format(experiment_name))
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



#experiment_list = ['Hitt_Expri_253_alpha_2.5_epsilon_1_delta_0.001']
#path = '/mnt/data/simulations/paper/'
#path = '/mnt/data/simulations/paper/'
path = './'
vtkoutput = False
matplotlibexport = False
netgenrender = False
scale_factor = 0.001
mesh_switch = 'rect_edgy'
make_linftyplot = True
# experiment_list = os.listdir('/mnt/data/simulations/paper/simulations/{}'.format(mesh_switch))
experiment_list = ['LIGHT_Expri_NOEV_refined_200_alpha_1_epsilon_1_delta_0.0003']

if netgenrender == True:
    import netgen.gui

for experiment_name in experiment_list:
    print(experiment_name)
    if mesh_switch == 'center':
        pathdirect = path+'simulations/center/{}/'.format(experiment_name)
    elif mesh_switch == 'edgy':
        pathdirect = path+'simulations/edgy/{}/'.format(experiment_name)
    elif mesh_switch == 'obst':
        pathdirect = path+'simulations/obst/{}/'.format(experiment_name)
    elif mesh_switch == 'trap':
        pathdirect = path+'simulations/trap/{}/'.format(experiment_name)
    elif mesh_switch == 'rect':
        pathdirect = path+'simulations/rect/{}/'.format(experiment_name)
    elif mesh_switch == 'circ_hd':
        pathdirect = path+'simulations/circ_hd/{}/'.format(experiment_name)
    elif mesh_switch == 'rect_edgy':
        pathdirect = path+'simulations/rect_edgy/{}/'.format(experiment_name)
    else:
        pathdirect = path+'simulations/none/{}/'.format(experiment_name)

    print(pathdirect)
    try:
        os.mkdir(pathdirect+'vtk')
    except Exception as FileExistsError:
        try:
            shutil.rmtree(pathdirect+'vtk')
        except OSError as e:
            print ("Error: %s - %s." % (e.filename, e.strerror))

        os.mkdir(pathdirect+'vtk')

    timelist = reorder_timesteps(experiment_name,pathdirect)
    print(timelist)
    paramdict = opencsvfile(experiment_name,pathdirect)
    T = float(paramdict['T'])
    dt = float(paramdict['dt'])
    alpha = float(paramdict['alpha'])
    delta = float(paramdict['delta'])
    meshsize = float(paramdict['meshsize'])
    geometry = paramdict['geometry']
    order = int(paramdict['order'])
    midpointc = (0,0)
    radius = 1
    length = 0
    if mesh_switch == 'edgy':
        mesh,V,u,v,gfuL2 = meshgeneration_n_spaces_rec_edgy(meshsize,order,midpointc,radius)
    elif mesh_switch == 'center':
        mesh,V,u,v,gfuL2 = meshgeneration_n_spaces_rec_center(meshsize,order,midpointc,radius)
    elif mesh_switch == 'obst':
        mesh,V,u,v,gfuL2 = meshgeneration_n_spaces_rec_obst(meshsize,order,midpointc,radius)
    elif mesh_switch == 'trap':
        mesh,V,u,v,gfuL2 = meshgeneration_n_spaces_rec_trap(meshsize,order,length)
    elif mesh_switch == 'rect':
        mesh,V,u,v,gfuL2 = meshgeneration_n_spaces_rec_rect(meshsize,order,length)
    elif mesh_switch == 'rect_edgy':
        mesh,V,u,v,gfuL2 = meshgeneration_n_spaces_rec_rect_edgy(meshsize,order,length)
    elif mesh_switch == 'circ_hd':
        mesh,V,u,v,gfuL2 = meshgeneration_n_spaces_circle_hd(meshsize,order,midpointc,radius)
    else:
        mesh,V,u,v,gfuL2 = meshgeneration_n_spaces_rec(meshsize,order,midpointc,radius)

    gfucalc = GridFunction(V,name="u_n")

    Linftymin = []
    Linftymax = []
    t = []
    k = 0

    for i in timelist:
        if i[0] == 0.0:
            gfucalc.Load(str(i[1]))
            gfuL2.Set(gfucalc.components[0])
            # pdb.set_trace()
            x = []
            y = []
            z = []
            trigs = []
            for v in mesh.vertices:
                p = v.point
                x.append(p[0])
                y.append(p[1])
                z.append(gfucalc.components[0].vec[v.nr])
            linfmax = max(z)
            linfmin = min(z)
            if matplotlibexport == True:
                for e in mesh.Elements():
                    trigs.append( [v.nr for v in e.vertices] )
            if netgenrender == True:
                Draw(gfucalc)
                visoptions.scalfunction="u_n:1"
                visoptions.vecfunction = "None"
                visoptions.scaledeform1 = scale_factor
                visoptions.deformation = 1
            Linftymin.append(linfmin)
            Linftymax.append(linfmax)
            t.append(i[0])
            if vtkoutput == True:
                vtk = VTKOutput(ma=mesh,coefs=[gfucalc[0],gfucalc[1]],names=['rho','c'],filename=pathdirect+'vtk/{}_{}'.format(experiment_name,k),subdivision=2)
                vtk.Do()
        else:
            if k % 10 == 0:
                print(Integrate(gfucalc,mesh))
                gfucalc.Load(str(i[1]))
                gfuL2.Set(gfucalc.components[0])
                x = []
                y = []
                z = []
                trigs = []
                for v in mesh.vertices:
                    p = v.point
                    x.append(p[0])
                    y.append(p[1])
                    z.append(gfucalc.components[0].vec[v.nr])

                if matplotlibexport == True:
                    for e in mesh.Elements():
                        trigs.append( [v.nr for v in e.vertices] )

                linfmin_l2 = min(gfuL2.vec)
                linfmax_l2 = max(gfuL2.vec)
                Linftymin.append(linfmin_l2)
                Linftymax.append(linfmax_l2)
                t.append(float(i[0]))
                if netgenrender == True:
                    Redraw(True)


                # store vertex coorinates in x/y and function value in z array


                if min(z) < 0 or linfmin_l2 < 0:
                    print('>>>>>>>>>>>>>>>>>><this expri: {}'.format(experiment_name))
                #        break

                if vtkoutput == True:
                    vtk = VTKOutput(ma=mesh,coefs=[gfucalc[0],gfucalc[1]],names=['rho','c'],filename=pathdirect+'vtk/{}_{}'.format(experiment_name,k),subdivision=2)
                    vtk.Do()
                    print('Timestep: {}'.format(i))
                if matplotlibexport == True:
                    fig = plt.figure()
                    ax = fig.gca(projection='3d')
                    ax.grid(True)
                    #fig.patch.set_facecolor('white')
                    surf = ax.plot_trisurf(x, y, z, triangles=trigs, antialiased=True, cmap=cm.coolwarm, linewidth=0, shade=False)
                    if alpha < 1:
                        zformatter = FormatStrFormatter('%.1f')
                        ax.w_zaxis.set_major_formatter(zformatter)
                    else:
                        pass

                    ax.set_xlabel('x')
                    ax.set_ylabel('y')
                    plt.savefig(path+'Experiments_Linftyplots/eps/{}_time_{}.eps'.format(experiment_name,str(i[0])))

                print(i[0])
        k += 1


    print(t)
    print(Linftymax)

    fig = plt.figure(1)
    plt.subplot(211)
    plt.plot(t,Linftymax)
    plt.title("Maximumvalues vs Time")
    plt.xlabel("t")
    plt.ylabel("Maximumvalue of rho")
    plt.subplot(212)
    plt.plot(t,Linftymin)
    plt.title("Minimumvalue vs Time")
    plt.xlabel("t")
    plt.ylabel("Minimumvalue of rho")
    # plt.plot(t[0:-2],Linftymax[0:-2])
    # plt.plot(t[0:-2],Linftymin[0:-2])
    #plt.draw()
    z = np.array([t,Linftymax])
    z= np.transpose(z)
    lasttime_plotted = round(T*dt,5)
    if make_linftyplot == True:
        df = pd.DataFrame(z, columns = ['t','f(t)'])
        df.to_csv(path+'Experiments_Linftyplots/{}_alpha_{}_delta_{}.dat'.format(experiment_name,alpha,delta), index=None, header=None)
        plt.savefig(path+'Experiments_Linftyplots/{}_endtime_{}.png'.format(experiment_name,i[0]))
    #plt.waitforbuttonpress(0)
    plt.close(fig)
