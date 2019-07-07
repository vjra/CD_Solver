from ngsolve import *
from netgen.geom2d import SplineGeometry
from ks_solver4_new_many_v2 import *
import csv
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import cm
import numpy as np
from sys import argv
from matplotlib.tri import Triangulation, TriAnalyzer, UniformTriRefiner


netgenswitch = False
timestep_to_plot = argv[1:]
experiment_name = 'LIGHT_Expri_NOEV_refined_200_alpha_1_epsilon_1_delta_0.0003'
# experiment_name = 'DIA_Expri_HD_102_alpha_1_epsilon_1_delta_0.005'
# experiment_name = 'DIA_Expri_HD_103_alpha_1_epsilon_1_delta_0.001'
# experiment_name = 'DIA_Expri_HD_503_alpha_1_epsilon_1_delta_0.0025'
# experiment_name = 'DIA_Expri_HD_104_alpha_1_epsilon_1_delta_0.0005'
mesh_switch = 'rect_edgy'
# 0 is cell density, concentration of chemoattractant is 1
density_concentration_switch = 1
# path = '/mnt/data/simulations/paper/'
path = './'
# path = '/mnt/data/simulations/tba/'
# path = '/mnt/data/simulations/paper/'
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

print(pathdirect)
for i in timestep_to_plot:
    print(i)
    print(pathdirect+'{}_time_{}'.format(experiment_name,i))
    gfucalc.Load(pathdirect+'{}_time_{}'.format(experiment_name,i))
    if netgenswitch == True:
        import netgen.gui
        Draw(gfucalc.components[density_concentration_switch],mesh,'density')
        visoptions.scalfunction="density"
        visoptions.vecfunction = "None"
        visoptions.scaledeform1 = 0.001
        visoptions.deformation = 0

    gfuL2.Set(gfucalc.components[density_concentration_switch])
    linfmin = min(gfuL2.vec)
    linfmax = max(gfuL2.vec)
    x = []
    y = []
    z = []
    trigs = []
    gfuL2.Set(gfucalc.components[density_concentration_switch])
    linfmin_l2 = min(gfuL2.vec)
    linfmax_l2 = max(gfuL2.vec)


    # store vertex coorinates in x/y and function value in z array
    for v in mesh.vertices:
        p = v.point
        x.append(p[0])
        y.append(p[1])
        z.append(gfucalc.components[density_concentration_switch].vec[v.nr])


    # triangulation
    for e in mesh.Elements():
        trigs.append( [v.nr for v in e.vertices] )

    # trigo = Triangulation(x, y, triangles=trigs)
    # refiner = UniformTriRefiner(trigo)
    # tri_refi, z_test_refi = refiner.refine_field(z, subdiv=2)

    fig, ax = plt.subplots()
    ax.grid(False)
    #fig.patch.set_facecolor('white')

    levels = np.linspace(0,max(z),10)
    levels = np.sort(np.append(levels,1))
    tcs = ax.tricontourf(x,y,trigs,z, levels=levels,cmap = cm.coolwarm)
    plt.ylim(0, 0.4)
    cbar = fig.colorbar(tcs)
    plt.xlabel('x')
    plt.ylabel('y')
    # fig.colorbar(surf, shrink=0.5, aspect=10)
    # surf.set_edgecolor(surf.to_rgba(surf._A)) # probably transparent edges
    # surf.set_edgecolor('black')
    # x = surf.get_edgecolors()
    # print(x)
    if density_concentration_switch == 0:
        plt.savefig('/mnt/data/simulations/DISS/paper/Experiments_Linftyplots/eps/HM_{}_time_{}.eps'.format(experiment_name,round(float(i),5)))
    elif density_concentration_switch == 1:
        plt.savefig('/mnt/data/simulations/DISS/paper/Experiments_Linftyplots/eps/HM_c_{}_time_{}.eps'.format(experiment_name,round(float(i),5)))
    else:
        print('Choose cell or chemoattractant density')

    print(min(z))
    print(max(z))
    print(linfmin_l2)
    print(linfmax_l2)
    print('done')
    plt.draw()
    plt.waitforbuttonpress(0)
    plt.close(fig)
