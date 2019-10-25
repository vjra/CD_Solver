from experiment_class import *
import netgen.gui

############################################################################
########################## Experiment class info ###########################
############################################################################

# Solves the parabolic-parabolic system for u=(\rho,c) \in \R^2
#
# u_t = div(A(u) \nabla u) + f(u)
#
# for t > 0 on a bounded domain Omega,
# with source term
# f(u) = (0,-u[1]+u[0]^alpha)
# and parameters
# with eps >= 0 and alpha >0.

############################################################################
########################## New features, not implemented ###################
############################################################################
# Solves the parabolic-parabolic system for u=(\rho,c) \in \R^2
#
# u_t = div(A(u) \nabla u) + f(u)
#
# for t > 0 on a bounded domain Omega,
# with source term
# f(u) = p(u)
# with p a polynomial and parameters
# with eps >= 0 and alpha >0.



filedir ='.'
# simulations_folder = '/simulations_data/'
simulations_folder = '/simulations_data/'
path = filedir+simulations_folder
visualoutput_solver = True
ini_data = (8*pi/(pi*1/200)*exp(-200*(x-0.5)**2-200*(y-0.5)**2), 0)
ini_data_str = '(8*pi/(pi*1/200)*exp(-200*(x-0.5)**2-200*(y-0.5)**2), 0)'
#expri_name,alpha,delta,epsilon,dt,T,meshsize,order,geometry,visualoutput_solver):
experiment1 = experiment('Test_version',1,10**(-1),0,10**(-2),0.05,0.03,2,'Unit square',visualoutput_solver)

experiment1.set_file_structure(simulations_folder,path,filedir)
experiment1.meshgeneration_n_spaces()
A = (1,-experiment1.u[0],experiment1.expri_data['delta'],1)
experiment1.set_initial_data(ini_data,ini_data_str)
experiment1.weak_formulation(A)
experiment1.plotseries_and_more(make_linftyplot = True,deform_param = 1, netgenrender = True)
# experiment1.plot_and_export([0.0,0.03,0.04],0,netgenrender = True,matplotlib_export = True)
