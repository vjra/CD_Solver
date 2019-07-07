from ngsolve import *
from netgen.geom2d import SplineGeometry
from ks_solver5_v1 import *
import csv
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import cm
import numpy as np
from sys import argv
import matplotlib.tri as tri


netgenswitch = True
timestep_to_plot = argv[1:]
experiment_name = 'DIA_Expri_HD_101_alpha_1_epsilon_1_delta_0.01'
mesh_switch = 'circ'
# path = '/mnt/data/simulations/paper/'
path = './'
# path = '/mnt/data/simulations/tba/'
#path = '/mnt/data/simulations/paper/'
def opencsvfile(name,pathdirect):
    with open(pathdirect+'parameter.csv', 'r') as csvfile:
        csv_reader = csv.DictReader(csvfile)
        paramdict = list(csv_reader)

    return paramdict[0]


if mesh_switch == 'center':
    pathdirect = path+'simulations/center/{}/'.format(experiment_name)
elif mesh_switch == 'edgy':
    pathdirect = path+'simulations/edgy/{}/'.format(experiment_name)
elif mesh_switch == 'obst':
    pathdirect = path+'simulations/obst/{}/'.format(experiment_name)
elif mesh_switch == 'circ':
    pathdirect = path+'simulations/circ_hd/{}/'.format(experiment_name)
else:
    pathdirect = path+'simulations/none/{}/'.format(experiment_name)


paramdict = opencsvfile(experiment_name,pathdirect)
T = float(paramdict['T'])
dt = float(paramdict['dt'])
alpha = float(paramdict['alpha'])
delta = float(paramdict['delta'])
meshsize = float(paramdict['meshsize'])
geometry = paramdict['geometry']
order = int(paramdict['order'])
midpointc = (0,0)
############################ATTENTION RADIUS###################################
radius = 1
if mesh_switch == 'center':
    mesh,V,u,v,gfuL2 = meshgeneration_n_spaces_rec_center(meshsize,order,midpointc,radius)
elif mesh_switch == 'edgy':
    mesh,V,u,v,gfuL2 = meshgeneration_n_spaces_rec_edgy(meshsize,order,midpointc,radius)
elif mesh_switch == 'obst':
    mesh,V,u,v,gfuL2 = meshgeneration_n_spaces_rec_obst(meshsize,order,midpointc,radius)
elif mesh_switch == 'circ':
    mesh,V,u,v,gfuL2 = meshgeneration_n_spaces_circle_hd_c(meshsize,order,midpointc,radius)
else:
    mesh,V,u,v,gfuL2 = meshgeneration_n_spaces_rec(meshsize,order,midpointc,radius)

gfucalc = GridFunction(V,name="u_n")
plt.ion()
print(pathdirect)
for i in timestep_to_plot:
    print(i)
    print(pathdirect+'{}_time_{}'.format(experiment_name,i))
    gfucalc.Load(pathdirect+'{}_time_{}'.format(experiment_name,i))
    if netgenswitch == True:
        import netgen.gui
        Draw(gfucalc.components[0])
        # visoptions.scalfunction="u_n:1"
        # visoptions.vecfunction = "None"
        # visoptions.scaledeform1 = 0.1
        # visoptions.deformation = 1

    gfuL2.Set(gfucalc.components[0])
    linfmin = min(gfuL2.vec)
    linfmax = max(gfuL2.vec)
    x = []
    y = []
    z = []
    trigs = []
    gfuL2.Set(gfucalc.components[0])
    linfmin_l2 = min(gfuL2.vec)
    linfmax_l2 = max(gfuL2.vec)


    # store vertex coorinates in x/y and function value in z array
    for v in mesh.vertices:
        p = v.point
        x.append(p[0])
        y.append(p[1])
        z.append(gfucalc.components[0].vec[v.nr])


    # triangulation
    for e in mesh.Elements():
        trigs.append( [v.nr for v in e.vertices] )

    cmap = cm.get_cmap(name='Blues', lut=None)
    trigo = tri.Triangulation(x, y, triangles=trigs)
    print('hi')
    fig, ax = plt.subplots()
    # ax.triplot(trigo, lw=0.5, color='white')
    # ax = fig.gca(projection='3d')
    # tcf = ax.tricontourf(trigs, z)
    # fig.colorbar(tcf)
    ax.grid(False)
    #fig.patch.set_facecolor('white')
    # surf = ax.plot_trisurf(x, y, z, triangles=trigs, antialiased=False, cmap=cm.coolwarm, linewidth=0, shade=False)
    # levels = np.arange(0,max(z), 2)
    # ax.tricontourf(x,y,trigs, z,levels=levels)
    # ax.tricontourf(x,y,trigs, z)
    tcs = ax.tricontour(x,y,trigs,z,levels = [10**(-2)], cmap = None, linewidths=1)
    ax.clabel(tcs, fontsize=10)
    # plt.savefig('./Experiments_Linftyplots/eps/{}_time_{}.eps'.format(experiment_name,round(float(i),5)))
    print(min(z))
    print(max(z))
    print(linfmin_l2)
    print(linfmax_l2)
    print('done')
    plt.show()
    # plt.waitforbuttonpress(0)
    # plt.close(fig)
