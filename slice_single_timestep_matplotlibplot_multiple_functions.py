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


netgenswitch = False
timestep_to_plot = argv[1:]
experiment_list = os.listdir('/home/oleingan/Nextcloud/KS2/ngsolve/repos/KS_ngsolve/simulations/circ_hd')
# experiment_list = ['DIA_Expri_HD_101_alpha_1_epsilon_1_delta_0.01','DIA_Expri_HD_102_alpha_1_epsilon_1_delta_0.005','DIA_Expri_HD_103_alpha_1_epsilon_1_delta_0.001','DIA_Expri_HD_104_alpha_1_epsilon_1_delta_0.0005']
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






############################ATTENTION RADIUS###################################



plt.ion()
fig, ax = plt.subplots()
ax.grid(False)
circle_boundary = plt.Circle((0, 0), 1, fill=False)
ax.add_artist(circle_boundary)
for expri in experiment_list:
    experiment_name = expri
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
    print(dt)
    if dt == 0.01:
        i = 4.99
    else:
        i = 4.999
    alpha = float(paramdict['alpha'])
    delta = float(paramdict['delta'])
    meshsize = float(paramdict['meshsize'])
    geometry = paramdict['geometry']
    order = int(paramdict['order'])
    midpointc = (0,0)
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
    print(pathdirect)
    gfucalc = GridFunction(V,name="u_n")
    print(i)
    print(pathdirect+'{}_time_{}'.format(experiment_name,i))
    gfucalc.Load(pathdirect+'{}_time_{}'.format(experiment_name,i))

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

    # cmap = cm.get_cmap(name='Blues', lut=None)
    trigo = tri.Triangulation(x, y, triangles=trigs)

    tcs = ax.tricontour(x,y,trigs,z,levels = [10**(-2)],colors = 'b', cmap = None, linewidths=1)
    fmt = {}
    deltastring = str(delta)
    fmt[10**(-2)] = r'$\delta=$'+deltastring
    print(fmt)
    ax.clabel(tcs,fmt=fmt, inline=1, fontsize=8)
    # ax.scatter(x,y,c = z)
    xnump = np.array(x)
    lineo = xnump
    ax.plot(x,lineo)
    print(min(z))
    print(max(z))
    print(linfmin_l2)
    print(linfmax_l2)
    print('done')
    input('what for next one')
    plt.pause(0.05)
    plt.show()

plt.savefig('./Experiments_Linftyplots/eps/delta_rad.eps')
    # plt.waitforbuttonpress(0)
    # plt.close(fig)
