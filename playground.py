import netgen.gui
from ngsolve import *
from math import pi
from ngsolve.internal import visoptions
from time import time
from ks_solver4_new_many import *
from netgen.geom2d import SplineGeometry

midpointr = (0,0)
radius = 1
midpointc = (0,0)
meshsize = 0.035
geometry = 'unitcircle'
order = 2
# Bumpfunction parameters
# Bumpfunction parameters
    # Bumpfunction parameters
mass1 = 5*pi
mass2 = 2*pi
mass3 = 2*pi
xi = 1/200
xi2 = 0

mesh,V,u,v,gfuL2 = meshgeneration_n_spaces_rec_edgy(meshsize,order,midpointc,radius)
gfucalc = initial_data_many(mass1,V,midpointr,xi,mass2,mass3)

Draw(gfucalc)
visoptions.vecfunction = "None"
visoptions.scaledeform1 = 0.01
visoptions.deformation = 1
Redraw(True)
