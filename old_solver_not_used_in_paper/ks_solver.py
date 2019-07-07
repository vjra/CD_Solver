from netgen.geom2d import unit_square
from ngsolve import *
from math import pi
from ngsolve.internal import visoptions
from time import time
from netgen.geom2d import SplineGeometry
import os
import csv
import numpy as np

# Solve (u_t,eps*v_t) = div(A*(u,v)) + (0,-u[1]+u[0]^alpha)

# filedir name where the simulations have to be saved
filedir = './simulations/'
visualoutput_solver = True
def meshgeneration_n_spaces(meshsize,polyorder):
    mesh = Mesh (unit_square.GenerateMesh(maxh=meshsize))
    V = VectorH1(mesh, order=polyorder, dirichlet=[])
    u,v = V.TnT()
    l2 = L2(mesh, order=0)
    gfuL2 = GridFunction(l2)
    return mesh,V,u,v,gfuL2

def meshgeneration_n_spaces_rec(meshsize,polyorder,midpoint,radius):
    geo = SplineGeometry()
    geo.AddCircle(midpoint, radius)
    mesh = Mesh (geo.GenerateMesh(maxh=meshsize))
    V = VectorH1(mesh, order=polyorder, dirichlet=[])
    u,v = V.TnT()
    l2 = L2(mesh, order=0)
    gfuL2 = GridFunction(l2)
    return mesh,V,u,v,gfuL2

def weak_formulation(uvold,V,u,v,dt,delta,epsilon,alpha):
    A = CoefficientFunction((1,-u[0],delta,1),dims=(2,2))
    a = BilinearForm(V)
    a += SymbolicBFI( u[0]*v[0] + epsilon*u[1]*v[1]+dt*InnerProduct(A*grad(u),grad(v)) + dt*CoefficientFunction((0,u[1]-u[0].Norm()**alpha))*v - uvold[0] * v[0]-epsilon*uvold[1]*v[1])
    return a

def initial_data(V,xi,xi2,mass,mass2,midpointr,midpointc):
    uvold = GridFunction(V)
    uvold.Set( CoefficientFunction((mass/(pi*xi)*exp(-1/xi*(x-midpointr[0])**2-1/xi*(y-midpointr[1])**2), (mass2/(pi*xi2)*exp(-1/xi2*(x-midpointc[0])**2-1/xi2*(y-midpointc[1])**2)  ))))
    return uvold

def initial_data_radial(V,xi,xi2,mass,mass2,midpointr,midpointc):
    uvold = GridFunction(V)
    uvold.Set( CoefficientFunction((mass/(pi*xi)*exp(-1/xi*(x-midpointr[0])**2-1/xi*(y-midpointr[1])**2),0)))
    return uvold

def initial_data_hitt(V):
    uvold = GridFunction(V)
    uvold.Set(CoefficientFunction(((80*(x**2+y**2-1)**2*(x-0.1)**2+5),0)))
    return uvold

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


def paramlist(experiment_name,xi,xi2,mass,mass2,alpha,epsilon,delta,dt,meshsize,order,T,geometry,midpointr,midpointc):
    paramdict = {'experiment_name': experiment_name,'xi': xi,'xi2': xi2,'mass': mass,'mass2': mass2,'alpha': alpha, 'epsilon': epsilon, 'delta': delta, 'dt': dt,'meshsize': meshsize,'order': order, 'T': T, 'geometry': geometry,'midpointr': midpointr, 'midpointc':midpointc}
    try:
        os.mkdir(filedir+experiment_name)
    except Exception as FileExistsError:
        pass

    with open(filedir+'{}/parameter.csv'.format(experiment_name),'w') as csvfile:
        fieldnames = ['experiment_name','xi','xi2','mass','mass2','alpha','epsilon','delta','dt','meshsize','order','T','geometry','midpointr','midpointc']
        paramwriter = csv.DictWriter(csvfile, fieldnames = fieldnames)
        paramwriter.writeheader()
        paramwriter.writerow(paramdict)

    return paramdict


def run(experiment_name,uvold,gfucalc,gfuL2,a,mesh,dt,nsteps,paramlisto,startingtimestep):
    try:
        os.mkdir(filedir+experiment_name)
    except Exception as FileExistsError:
        pass

    initialmass = Integrate(uvold,mesh)
    with TaskManager():
        for i in range(startingtimestep,nsteps):
            atime = time()
            uvold.Set(gfucalc)
            SimpleNewtonSolve(gfucalc,a)
            mass = Integrate(gfucalc,mesh)
            gfucalc.Save(filedir+'{}/{}_timestep_{}'.format(experiment_name,experiment_name,i))
            with open(filedir+'{}/last_time_step.txt'.format(experiment_name),"w") as f:
                f.write(str(int(i)))

            gfuL2.Set(gfucalc.components[0])
            linfmin = min(gfuL2.vec)
            linfmax = max(gfuL2.vec)
            print("total mass = ", mass,"Minvalue: ", linfmin, "Maxvalue :", linfmax)
            if mass[0] > initialmass[0]*2 or linfmin < -10:
                errormess(mass)
                with open(filedir+'{}/errormessage.txt'.format(experiment_name),"w") as f:
                    f.write('''Error at timestep {}:\nMass {}
                               \n min: {} \n max: {}'''.format(i,mass,linfmin,linfmax))

                with open(filedir+'{}/error_endtime.txt'.format(experiment_name),"w") as f:
                    f.write(str(int(i)))

                break
            if visualoutput_solver == True:
                Redraw(True)
            else:
                pass
            btime = time()
            print("step took: ",str(btime-atime))
            print("Time: {}, time step: {} of {}".format(str(dt*i),i,nsteps))

        with open(filedir+'{}/parameter.txt'.format(experiment_name),"w") as f:
            for param in paramlisto:
                f.write(str(param) + ': '+str(paramlisto[param])+'\n')
    return gfucalc
