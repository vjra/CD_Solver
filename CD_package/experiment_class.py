from solver import *
from ngsolve import *
import netgen.gui
# Solve (u_t,eps*v_t) = div(A*(u,v)) + (0,-u[1]+u[0]^alpha)

############################### global variables #################################
# subfolder where the simulations will be stored
filedir ='.'
simulations_folder = '/simulations_data/'
path = filedir+simulations_folder



class experiment():
    def __init__(self,expri_name,alpha,delta,epsilon,dt,T,meshsize,order,geometry,visualoutput_solver):
        self.expri_data = {'expri_name': expri_name,'alpha': alpha, 'delta': delta, \
                           'epsilon': epsilon, 'dt': dt,'T': T, 'meshsize': meshsize,'order':order,'geometry': geometry}
        self.experiment_full_name = expri_name+'_{}_alpha_{}_epsilon_{}_delta_{}'.format(geometry.replace(" ", ""),alpha,epsilon,delta)
        self.visualoutput_solver = visualoutput_solver
        self.simulations_folder = simulations_folder
        self.path = path
        self.filedir = filedir

    def set_file_structure(self,simulations_folder,path,filedir):
        self.simulations_folder = simulations_folder
        self.path = path
        self.filedir = filedir

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

    def set_initial_data(self,ini_data,ini_data_str='No string was entered'):
        self.uvold = GridFunction(self.FEspace)
        self.ini_data_str = '{}'.format(ini_data_str)
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

    def run_experiment(self):
        experiment_full_name = self.experiment_full_name
        mesh = self.mesh
        FEspace = self.FEspace
        weakform = self.a
        u = self.u
        v = self.v
        gfuL2 = self.gfuL2
        gfucalc = self.gfucalc
        dt = self.expri_data['dt']
        alpha = self.expri_data['alpha']
        delta = self.expri_data['delta']
        epsilon = self.expri_data['epsilon']
        T = self.expri_data['T']
        ini_data_str = self.ini_data_str
        meshsize = self.expri_data['meshsize']
        order = self.expri_data['order']
        geometry = self.expri_data['geometry']
        visualoutput_solver = self.visualoutput_solver
        if visualoutput_solver == True:
            import netgen.gui

        simulations_folder = self.simulations_folder
        path = self.path
        filedir = self.filedir

        print('simulations saved at: '+path+experiment_full_name)
        try:
            os.mkdir(path+experiment_full_name)
        except Exception as FileExistsError:
            pass

        cont_switch = folder_checker(experiment_full_name,simulations_folder)
        print('Contswitch variable: '+str(cont_switch))
        if cont_switch == True:
            old_timestepping_last_time, uvold = continuation_time(experiment_full_name,path,FEspace,dt)
            startingtimestep = int(old_timestepping_last_time/dt)+1
        elif cont_switch == False:
            uvold = self.uvold
            uvold.Save(filedir+'{}/{}_time_{}'.format(experiment_full_name,experiment_full_name,0))
            startingtimestep = 1
        else:
            print('something went wrong in cont_switching')
        mass = (Integrate(uvold.components[0],mesh),Integrate(uvold.components[1],mesh))
        gfuL2.Set(uvold.components[0])
        linfmin_l2 = min(gfuL2.vec)
        linfmax_l2 = max(gfuL2.vec)
        nsteps = int(floor(T/dt))
        print('Number of time steps: ', nsteps)
        print('Total mass of initial data: ', mass)
        if visualoutput_solver == True:
            gfucalc.vec.data = uvold.vec
            Draw(gfucalc.components[1],mesh, 'c')
            Draw(gfucalc.components[0],mesh, name = 'rho')
            # visoptions.scalfunction="rho"
            # visoptions.vecfunction = "None"
            # visoptions.scaledeform1 = 0.001
            # visoptions.deformation = 1
            # Redraw(True)
        else:
            gfucalc.vec.data = uvold.vec
        time_list_dict = {'dt': dt, 'timestep': 0,'filename':simulations_folder+experiment_full_name+'/'+experiment_full_name+'_time_0',\
                          'linfymin': round(linfmin_l2,3),'linfmax':  round(linfmax_l2,3),'mass1': round(mass[0],3),'mass2':round(mass[1],3)}
        with open(path+'{}/time_list.csv'.format(experiment_full_name),'w') as csvfile:
            fieldnames = ['dt','timestep','filename','linfymin','linfmax','mass1','mass2']
            listwriter = csv.DictWriter(csvfile, fieldnames = fieldnames)
            listwriter.writerow(time_list_dict)


        sleep(1)
        log_file_name = log_file_creator(simulations_folder,experiment_full_name)
        paramlisto = paramlist(path,experiment_full_name,mass,alpha,epsilon,delta,dt,meshsize,order,T,geometry,ini_data_str)
        SetNumThreads(8)
        sleep(1)
        self.gfucalc,endtime = run(path,experiment_full_name,uvold,gfucalc,gfuL2,weakform,mesh,dt,nsteps,paramlisto,startingtimestep,visualoutput_solver)
        print('Experiment: {}, done'.format(experiment_full_name))

    def plot_and_more(self):
