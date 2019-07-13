from netgen.geom2d import unit_square
from ngsolve import *
from math import pi
from ngsolve.internal import visoptions
from time import time
from netgen.geom2d import SplineGeometry
import os
import csv
import numpy as np
from netgen.meshing import PointId
import netgen.gui
# Bumpfunction parameters
mass1 = 0.05
n_of_drops = 2000
xi = 1/800
length = 0
order = 2
# number_of_mesh_elements_for_refinement
# thickness_refinement_layer

class experiment():
    def __init__(self,expri_name,alpha,delta,epsilon,dt,T,meshsize,order,geometry):
        self.expri_data = {'expri_name': expri_name,'alpha': alpha, 'delta': delta, \
                           'epsilon': epsilon, 'dt': dt,'T': T, 'meshsize': meshsize,'order':order,'geometry': geometry}

        self.experiment_full_name = expri_name+'_{}_alpha_{}_epsilon_{}_delta_{}'.format(geometry.replace(" ", ""),alpha,epsilon,delta)

    def meshgeneration_n_spaces(self):
        meshsize = self.expri_data['meshsize']
        order = self.expri_data['order']
        self.mesh = Mesh(unit_square.GenerateMesh(maxh=meshsize))
        X = H1(self.mesh, order=order)
        Q = H1(self.mesh, order=order)
        self.FEspace = FESpace([X,Q])
        self.u = self.FEspace.TrialFunction()
        self.v = self.FEspace.TestFunction()
        l2 = L2(self.mesh, order=0)
        self.gfuL2 = GridFunction(l2)
        self.gfucalc = GridFunction(self.FEspace, name = "u_h")

    def set_initial_data(self,ini_data):
        self.uvold = GridFunction(self.FEspace)
        self.uvold.components[0].Set(CoefficientFunction(ini_data[0]))
        self.uvold.components[1].Set(CoefficientFunction(ini_data[1]))

    def weak_formulation(self,diffusion_matrix,source_term=0):
        self.diffusion_matrix = CoefficientFunction(diffusion_matrix,dims=(2,2))
        self.a = BilinearForm(self.FEspace)
        u = self.u
        v = self.v
        dt = self.expri_data['dt']
        alpha = self.expri_data['alpha']
        delta = self.expri_data['delta']
        epsilon = self.expri_data['epsilon']
        uvold = self.uvold
        gradu = CoefficientFunction( ( grad(self.u[0]), grad(self.u[1]) ), dims = (2,2) )
        gradv = CoefficientFunction( ( grad(v[0]), grad(v[1]) ), dims = (2,2) )
        self.a += SymbolicBFI(u[0]*v[0] + epsilon*u[1]*v[1]+dt*InnerProduct(self.diffusion_matrix*gradu,gradv) + dt*CoefficientFunction( (0,u[1]-u[0].Norm()**alpha) )*CoefficientFunction( (v[0],v[1]) ) - uvold.components[0] * v[0]-epsilon*uvold.components[1]*v[1])


ini_data = (8*pi/(pi*1/200)*exp(-200*(x-0.5)**2-200*(y-0.5)**2), 0)
experiment1 = experiment('Test_version',1,0,0,10**(-2),1,0.04,2,'Unit square')
# A = (1,-experiment1.u[0],experiment1.delta,1)
A = (1,0,0,1)
experiment1.meshgeneration_n_spaces()
experiment1.set_initial_data(ini_data)
print(experiment1.u)
experiment1.weak_formulation(A)
print(experiment1.diffusion_matrix)
# Draw(experiment1.uvold)
# Draw(experiment1.uvold.components[0])
# # visoptions.scalfunction="u_n:1"
# visoptions.vecfunction = "None"
# visoptions.scaledeform1 = 0.0005
# visoptions.deformation = 1

#

#

#
#
#
#
#
# def SimpleNewtonSolve(gfu,a,tol=1e-13,maxits=25):
#     res = gfu.vec.CreateVector()
#     du = gfu.vec.CreateVector()
#     fes = gfu.space
#     for it in range(maxits):
#         print ("Iteration {:3}  ".format(it),end="")
#         a.Apply(gfu.vec, res)
#         a.AssembleLinearization(gfu.vec)
#         du.data = a.mat.Inverse(fes.FreeDofs()) * res
#         gfu.vec.data -= du
#
#         #stopping criteria
#         stopcritval = sqrt(abs(InnerProduct(du,res)))
#         print ("<A u",it,", A u",it,">_{-1}^0.5 = ", stopcritval)
#         if stopcritval < tol:
#             break
#
# def errormess(mass):
#     print("ERROR: Minimumvalue exploading or mass negative")
#     print("Mass: {}".format(mass))
#
#
# def paramlist(experiment_name,xi,xi2,mass,mass2,alpha,epsilon,delta,dt,meshsize,order,T,geometry,midpointr,midpointc):
#     paramdict = {'experiment_name': experiment_name,'xi': xi,'xi2': xi2,'mass': mass,'mass2': mass2,'alpha': alpha, 'epsilon': epsilon, 'delta': delta, 'dt': dt,'meshsize': meshsize,'order': order, 'T': T, 'geometry': geometry,'midpointr': midpointr, 'midpointc':midpointc}
#     try:
#         os.mkdir(filedir+experiment_name)
#     except Exception as FileExistsError:
#         pass
#
#     with open(filedir+'{}/parameter.csv'.format(experiment_name),'w') as csvfile:
#         fieldnames = ['experiment_name','xi','xi2','mass','mass2','alpha','epsilon','delta','dt','meshsize','order','T','geometry','midpointr','midpointc']
#         paramwriter = csv.DictWriter(csvfile, fieldnames = fieldnames)
#         paramwriter.writeheader()
#         paramwriter.writerow(paramdict)
#
#     return paramdict
#
#
# def time_list_writer(experiment_name,dt,timestep,filename,linfmin_l2,linfmax_l2,mass):
#     time_list_dict = {'dt': dt, 'timestep': timestep,'filename':filename,'linfymin': round(linfmin_l2,3),'linfmax':  round(linfmax_l2,3),'mass1': round(mass[0],3),'mass2':round(mass[1],3)}
#
#     with open(filedir+'{}/time_list.csv'.format(experiment_name),'a') as csvfile:
#         fieldnames = ['dt','timestep','filename','linfymin','linfmax','mass1','mass2']
#         listwriter = csv.DictWriter(csvfile,fieldnames = fieldnames)
#         listwriter.writerow(time_list_dict)
#
#
#
# def run(experiment_name,uvold,gfucalc,gfuL2,a,mesh,dt,nsteps,paramlisto,startingtimestep,visualoutput_solver):
#
#     initialmass = (Integrate(uvold.components[0],mesh),Integrate(uvold.components[1],mesh))
#     with TaskManager():
#         for i in range(startingtimestep,nsteps):
#             atime = time()
#             uvold.vec.data = gfucalc.vec
#             # uvold.Set(gfucalc)
#             print('inside')
#             SimpleNewtonSolve(gfucalc,a)
#             mass = (Integrate(gfucalc.components[0],mesh),Integrate(gfucalc.components[1],mesh))
#             gfuL2.Set(gfucalc.components[0])
#             linfmin_l2 = min(gfuL2.vec)
#             linfmax_l2 = max(gfuL2.vec)
#             modulo_constant = 1
#             if dt < 1e-5:
#                 modulo_constant = 100
#             elif dt < 1e-9:
#                 modulo_constant = 1000
#
#             if i % modulo_constant == 0:
#                 gfucalc.Save(filedir+'{}/{}_time_{}'.format(experiment_name,experiment_name,str(i*dt)))
#                 with open(filedir+'{}/last_time.txt'.format(experiment_name),"w") as f:
#                     f.write(str(i*dt))
#                 with open(filedir+'{}/last_time_before.txt'.format(experiment_name),"w") as f:
#                     f.write(str((i-modulo_constant)*dt))
#                 filename_saved = filedir+'{}/{}_time_{}'.format(experiment_name,experiment_name,str(i*dt))
#                 print(experiment_name,dt,i,filename_saved)
#                 time_list_writer(experiment_name,dt,i,filename_saved,linfmin_l2,linfmax_l2,mass)
#
#
#             if  linfmin_l2 < 0:
#                 print('>>>>>>>>>>>>>>>>> Negative Value >>>>>>>>>>>>>>>>> {}'.format(linfmin_l2))
#                 break
#
#             print("total mass = ", mass, "MinvalueL2: ", linfmin_l2, "MaxvalueL2: ", linfmax_l2)
#             if mass[0] > initialmass[0]*2  or linfmin_l2 <0:
#                 endtime = i*dt
#                 errormess(mass)
#                 with open(filedir+'{}/errormessage.txt'.format(experiment_name),"w") as f:
#                     f.write('''Error at timestep {}:\nMass {}
#                                \n min: {} \n max: {}'''.format(i,mass,linfmin_l2,linfmax_l2))
#
#                 with open(filedir+'{}/error_endtime.txt'.format(experiment_name),"w") as f:
#                     f.write(str(int(i)))
#
#                 break
#             if visualoutput_solver == True:
#                 Redraw(True)
#
#             btime = time()
#             print("step took: ",str(btime-atime))
#             print("Time: {}, time step: {} of {}".format(str(dt*i),i,nsteps))
#             endtime = nsteps*dt
#
#         with open(filedir+'{}/parameter.txt'.format(experiment_name),"w") as f:
#             for param in paramlisto:
#                 f.write(str(param) + ': '+str(paramlisto[param])+'\n')
#
#
#     return gfucalc, endtime
