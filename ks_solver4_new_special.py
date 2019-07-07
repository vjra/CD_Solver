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

# Solve (u_t,eps*v_t) = div(A*(u,v)) + (0,-u[1]+u[0]^alpha)

# filedir name where the simulations have to be saved
filedir = './simulations/'
visualoutput_solver = True
# simulations paper:
# number_of_mesh_elements_edgy = 8000
# number_of_mesh_elements_center = 2800
# refinement_radius = 0.85
# refinement_radius_center = 0.4
number_of_mesh_elements_edgy = 12000
number_of_mesh_elements_center = 2800
refinement_radius = 0.85
refinement_radius_center = 0.4


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

def meshgeneration_n_spaces_rec_edgy(meshsize,polyorder,midpoint,radius):
    geo = SplineGeometry()
    geo.AddCircle(midpoint, radius)
    mesh = Mesh (geo.GenerateMesh(maxh=meshsize))
    bl = refinement_radius
    def Mark():
        for el in mesh.Elements():
            mesh.SetRefinementFlag(el, False)
        for el in mesh.Elements(BND):
            mesh.SetRefinementFlag(el, False)
        for el in mesh.Elements():
            for v in el.vertices:
                pnt = mesh.ngmesh[PointId(v.nr+1)]
                if (pnt[0]**2+pnt[1]**2>=bl):
                    mesh.SetRefinementFlag(el, True)

    while mesh.ne < number_of_mesh_elements_edgy:
        print(mesh.ne)
        Mark()
        mesh.Refine()
        print(mesh.nv, mesh.ne)

    V = VectorH1(mesh, order=polyorder, dirichlet=[])
    u,v = V.TnT()
    l2 = L2(mesh, order=0)
    gfuL2 = GridFunction(l2)
    return mesh,V,u,v,gfuL2

def meshgeneration_n_spaces_rec_center(meshsize,polyorder,midpoint,radius):
    geo = SplineGeometry()
    geo.AddCircle(midpoint, radius)
    mesh = Mesh (geo.GenerateMesh(maxh=meshsize))
    bl = refinement_radius_center
    def Mark():
        for el in mesh.Elements():
            mesh.SetRefinementFlag(el, False)
        for el in mesh.Elements(BND):
            mesh.SetRefinementFlag(el, False)
        for el in mesh.Elements():
            for v in el.vertices:
                pnt = mesh.ngmesh[PointId(v.nr+1)]
                if (pnt[0]**2+pnt[1]**2<=bl):
                    mesh.SetRefinementFlag(el, True)

    while mesh.ne < number_of_mesh_elements_center:
        print(mesh.ne)
        Mark()
        mesh.Refine()
        print(mesh.nv, mesh.ne)

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

def initial_data_many(mass1,V,midpointr,xi,mass2=0,mass3=0):
    uvold = GridFunction(V)
    uvold.Set(
              CoefficientFunction(
                                  (mass1/(pi*xi)*exp(-1/xi*(x-midpointr[0]-0.25)**2-1/xi*(y-midpointr[1])**2)
                                   +mass2/(pi*xi)*exp(-1/xi*(x-midpointr[0]+0.25)**2-1/xi*(y-midpointr[1])**2)
                                   +mass3/(pi*xi)*exp(-1/xi*(x-midpointr[0])**2-1/xi*(y-midpointr[1]+0.25)**2)
                                   +mass3/(pi*xi)*exp(-1/xi*(x-midpointr[0])**2-1/xi*(y-midpointr[1]-0.25)**2)
                                   +mass3/(pi*xi)*exp(-1/xi*(x-midpointr[0])**2-1/xi*(y-midpointr[1]-0.5)**2)
                                   +mass3/(pi*xi)*exp(-1/xi*(x-midpointr[0])**2-1/xi*(y-midpointr[1]-0.35)**2)
                                   +mass3/(pi*xi)*exp(-1/xi*(x-midpointr[0]+0.5)**2-1/xi*(y-midpointr[1])**2)
                                   +mass3/(pi*xi)*exp(-1/xi*(x-midpointr[0]-0.5)**2-1/xi*(y-midpointr[1]-0.25)**2),0)))
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


def time_list_writer(experiment_name,dt,timestep,filename,linfmin_l2,linfmax_l2,mass):
    time_list_dict = {'dt': dt, 'timestep': timestep,'filename':filename,'linfymin': round(linfmin_l2,3),'linfmax':  round(linfmax_l2,3),'mass1': round(mass[0],3),'mass2':round(mass[1],3)}

    with open(filedir+'{}/time_list.csv'.format(experiment_name),'a') as csvfile:
        fieldnames = ['dt','timestep','filename','linfymin','linfmax','mass1','mass2']
        listwriter = csv.DictWriter(csvfile,fieldnames = fieldnames)
        listwriter.writerow(time_list_dict)



def run(experiment_name,uvold,gfucalc,gfuL2,a,mesh,dt,nsteps,paramlisto,startingtimestep,visualoutput_solver):

    initialmass = Integrate(uvold,mesh)
    with TaskManager():
        for i in range(startingtimestep,nsteps):
            atime = time()
            uvold.Set(gfucalc)
            SimpleNewtonSolve(gfucalc,a)
            mass = Integrate(gfucalc,mesh)
            gfuL2.Set(gfucalc.components[0])
            linfmin_l2 = min(gfuL2.vec)
            linfmax_l2 = max(gfuL2.vec)
            modulo_constant = 1
            if dt < 1e-5:
                modulo_constant = 100
            elif dt < 1e-9:
                modulo_constant = 1000

            if i % modulo_constant == 0:
                gfucalc.Save(filedir+'{}/{}_time_{}'.format(experiment_name,experiment_name,str(i*dt)))
                with open(filedir+'{}/last_time.txt'.format(experiment_name),"w") as f:
                    f.write(str(i*dt))
                with open(filedir+'{}/last_time_before.txt'.format(experiment_name),"w") as f:
                    f.write(str((i-modulo_constant)*dt))

                filename_saved = filedir+'{}/{}_time_{}'.format(experiment_name,experiment_name,str(i*dt))
                print(experiment_name,dt,i,filename_saved)
                time_list_writer(experiment_name,dt,i,filename_saved,linfmin_l2,linfmax_l2,mass)


            # if  linfmin_l2 < 0:
            #     print('>>>>>>>>>>>>>>>>> Negative Value >>>>>>>>>>>>>>>>> {}'.format(linfmin_l2))
            #
            # print("total mass = ", mass, "MinvalueL2: ", linfmin_l2, "MaxvalueL2: ", linfmax_l2)
            # if mass[0] > initialmass[0]*2  or linfmin_l2 <0:
            #     endtime = i*dt
            #     errormess(mass)
            #     with open(filedir+'{}/errormessage.txt'.format(experiment_name),"w") as f:
            #         f.write('''Error at timestep {}:\nMass {}
            #                    \n min: {} \n max: {}'''.format(i,mass,linfmin_l2,linfmax_l2))
            #
            #     with open(filedir+'{}/error_endtime.txt'.format(experiment_name),"w") as f:
            #         f.write(str(int(i)))
            #
            #     break
            if visualoutput_solver == True:
                Redraw(True)

            btime = time()
            print("step took: ",str(btime-atime))
            print("Time: {}, time step: {} of {}".format(str(dt*i),i,nsteps))
            endtime = i*dt

            # if dt == 1e-2 and dt*i > 1.5:
            #     print('jump')
            #     break

        with open(filedir+'{}/parameter.txt'.format(experiment_name),"w") as f:
            for param in paramlisto:
                f.write(str(param) + ': '+str(paramlisto[param])+'\n')




    return gfucalc, endtime