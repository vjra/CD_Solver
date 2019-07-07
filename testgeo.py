from netgen.geom2d import SplineGeometry
from ngsolve import *
from math import pi
from ngsolve.internal import visoptions
from ks_solver4_new_many_v2 import *
import netgen.gui


meshsize = 0.04
polyorder = 2
length = 0
xi = 1/800
xi2 = 0
mass = 20*pi
mass2 = 0
midpointc = (0,0)
midpointr = (0.5,0.2)

mesh,V,u,v,gfuL2 = meshgeneration_n_spaces_rec_rect(meshsize,polyorder,length)
uvold = initial_data_radial(V,xi,xi2,mass,mass2,midpointr,midpointc)
Draw(uvold)
