from experiment_class import *
import netgen.gui

############################################################################
##################### Example code to use CS_package #######################
############################################################################

# Solves parabolic-parabolic system for u=(\rho,c) \in \R^2
# u_t = div(A(u) \nabla u) + f(u)
#
# for t > 0 on a bounded domain Omega,
# with source term
# f(u) = (0,-u[1]+u[0]^alpha)
# and parameters
# with eps >= 0 and alpha >0.

############################################################################
################ Experiment initializing and file structure ################
############################################################################

# Define in which directory to save the simulation data.
filedir ='.'
# Define folder name for the simulation data.
simulations_folder = '/simulations_data/'
# Define path
path = filedir+simulations_folder
# Set if there should be a visual output via netgen.gui or not.
visualoutput_solver = True

# Define experiment object with parameters.

# experiment_name (str): used in all export or saving file.
# alpha, delta (float): Free parameters to
# epsilon (float): Parameter to switch from parabolic-parabolic to parabolic-elliptic modelself.
# dt (float): time step size.
# T (float): end time of the simulation.
# meshsize (float): mesh size.
# order (float): order of the functions in the finite element space.
# geometry (str): geometry specification as string.
# visualoutput_solver (boolean): True starts netgen for plotting during computations.

experiment1 = experiment('Test_version',1,10**(-2),0,10**(-2),0.5,0.03,2,'Unit square',visualoutput_solver)

experiment1.set_file_structure(simulations_folder,path,filedir)
# create spaces and mesh.
experiment1.meshgeneration_n_spaces()
# Define initial data as ngsolve CoefficientFunction.
ini_data = (8*pi/(pi*1/200)*exp(-200*(x-0.5)**2-200*(y-0.5)**2), 0)
# ini_data_str = '(8*pi/(pi*1/200)*exp(-200*(x-0.5)**2-200*(y-0.5)**2), 0)'

# define diffusion matrix, use experiment1.u[0] as proxy for \rho and
# experiment1.u[1] as proxy for c.
A = (1,-experiment1.u[0],experiment1.expri_data['delta'],1)
experiment1.set_initial_data(ini_data)
experiment1.weak_formulation(A)

############################################################################
############################## Solve experiment ############################
############################################################################
# run simulation
experiment1.run_experiment()
