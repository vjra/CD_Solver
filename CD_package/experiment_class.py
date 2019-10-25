from solver import *
from ngsolve import *
from plots_and_more_lib import *
import logging
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
########################## Global variables ################################
############################################################################


# subfolder where the simulations will be stored
filedir ='.'
simulations_folder = '/simulations_data/'
path = filedir+simulations_folder

############################################################################
########################## Logging #########################################
############################################################################

# create and configure logger
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

# create a file handler
handler = logging.FileHandler('./logs/experiment_class.log', mode='w')
handler.setLevel(logging.INFO)

# create a logging format
formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
handler.setFormatter(formatter)

# add the file handler to the logger
logger.addHandler(handler)

############################################################################
########################## Class definition ################################
############################################################################

class experiment():
    """

    """
    ############################################################################
    ########################## Class init and file structure ###################
    ############################################################################
    def __init__(self,experiment_name,alpha,delta,epsilon,dt,T,meshsize,order,geometry,visualoutput_solver):
        """ Initialize object of class experiment.
            In particular, it sets the name, parameters and filepath
            for the experiment.

        Attributes
        ----------

        experiment_name (str): used in all export or saving file.
        alpha, delta (float): Free parameters to
        epsilon (float): Parameter to switch from parabolic-parabolic to parabolic-elliptic modelself.
        dt (float): time step size.
        T (float): end time of the simulation.
        meshsize (float): mesh size.
        order (float): order of the functions in the finite element space.
        geometry (str): geometry specification as string.
        visualoutput_solver (boolean): True starts netgen for plotting during computations.

        Initilizing class variables
        ---------------------------

        expri_data (dict): Dictionary including simulation parameters and name.
            alpha, delta (float): Free parameters to
            epsilon (float): Parameter to switch from parabolic-parabolic to parabolic-elliptic modelself.
            dt (float): time step size.
            T (float): end time of the simulation.
            meshsize (float): mesh size.
            order (float): order of the functions in the finite element space.
            geometry (str): geometry specification as string.

        experiment_full_name (str): Full experiment name, including simulation parameter.
        visualoutput_solver (boolean): True starts netgen for plotting during computations.
        simulations_folder (str): Folder used to store simulation data.
        path (str): filedir/simulationfolder
        filedir (str): filedir of simulationfolder.

        Methods
        -------
        set_file_structure(self,simulations_folder,path,filedir):

        meshgeneration_n_spaces(self):

        set_initial_data(self,ini_data,ini_data_str='No string was entered'):

        weak_formulation(self,diffusion_matrix,source_term=0):

        run_experiment(self):

        plotseries_and_more(self, netgenrender = True, vtkexport = False,
                            make_linftyplot = False,scale_factor = 1,
                            plot_modulo = 1,deform_param = 1):

        plot_and_export(self, at_times = [0], sol_component = 0,netgenrender = True,
                            matplotlib_export = False, vtkexport = False, scale_factor = 1,
                            plot_modulo = 1,deform_param = 1):
        """
        logger.info("Creating experiment object: {}".format(experiment_name))
        self.expri_data = {'experiment_name': experiment_name,'alpha': alpha, 'delta': delta, \
                           'epsilon': epsilon, 'dt': dt,'T': T, 'meshsize': meshsize,'order':order,'geometry': geometry}
        self.experiment_full_name = experiment_name+'_{}_alpha_{}_epsilon_{}_delta_{}'.format(geometry.replace(" ", ""),alpha,epsilon,delta)
        self.visualoutput_solver = visualoutput_solver
        self.simulations_folder = simulations_folder
        self.path = path
        self.filedir = filedir

    def set_file_structure(self,simulations_folder,path,filedir):
        """ Can be used to change the file paths for a experiment object after
            creation.

        Parameters
        ----------

        simulations_folder (str): Folder used to store simulation data.
        path (str): filedir/simulationfolder
        filedir (str): filedir of simulationfolder.
        """
        logger.info("Setting file structure: simulation folder: {}, path: {}, filedir: {}"
                     .format(simulations_folder,path,filedir))
        self.simulations_folder = simulations_folder
        self.path = path
        self.filedir = filedir

    ############################################################################
    ############### Meshgeneration, Finite element spaces, #####################
    ###############    initial data and weak formulation   #####################
    ############################################################################
    def meshgeneration_n_spaces(self):
        """ Initializing finite element spaces for test and trial functions,
            using parameters of object.

        Initializing class variables
        ---------------------------

        FEspace (ngsolve.comp.FESpace object): Tensorproduct of X and Q, i.e.
                                               finite element space for \rho and c.
        u (list of ngsolve.comp.ProxyFunction): Symbolic function for weak formulation,
                                                u[0] symbolizes \rho and u[1] symbolizes c.
        v (list of ngsolve.comp.ProxyFunction): Symbolic function for weak formulation,
                                                v[0] symbolizes test function for first equation
                                                and v[1] symbolizes for second equation.
        gfucalc (ngsolve.comp.GridFunction): Gridfunction of order specified by attribute "order" and
                                             lifes in FEspace. Discrete function representing the solution of
                                             the computation. \rho = gfucalc.components[0], c =  gfucalc.components[1].
        gfuL2 (ngsolve.comp.GridFunction): Gridfunction of order 0. Using L2() from ngsolve. Used to approximate L^\infty norm.

        """
        logger.info("""Creating mesh and FE space with parameters: meshsize: {}, order: {}""".format(self.expri_data['meshsize'],self.expri_data['order']))
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
        """ Set initial data on gridfunction.

        Initializing class variables
        ---------------------------

        uvold (ngsolve.comp.GridFunction): Gridfunction of order specified by attribute "order" and
                                           lifes in FEspace. Discrete function representing the previous
                                           time step during computation.
                                           \rho_0 = uvold.components[0], c =  uvold.components[1].

        """
        logger.info("Setting Initial data.")
        self.uvold = GridFunction(self.FEspace)
        self.ini_data_str = '{}'.format(ini_data_str)
        self.uvold.components[0].Set(CoefficientFunction(ini_data[0]))
        self.uvold.components[1].Set(CoefficientFunction(ini_data[1]))

    def weak_formulation(self,diffusion_matrix,source_term=0):
        """ Defining weak formulation of the problem, using BilinearForm of ngsolve.

        Parameters
        ----------

        diffusion_matrix (tuple): Flattend diffusion matrix A for weak formulation.
                                  Length needs to be 4.
                                  A = a_11 a_12
                                      a_21 a_22
                                  => diffusion_matrix = (a_11,a_12,a_21,a_22)
                                  a_ij needs to be proxy function or CoefficientFunction object as in ngsolve.
                                  e.g.:
                                  1) A = (1,0,0,1) for Laplacian
                                  2) A = (1,-experiment1.u[0],0,1) for classical Keller--Segel
        source_term (tuple): Not yet implemented...


        Initializing class variables
        ---------------------------
        diffusion_matrix (ngsolve.fem.CoefficientFunction): Defining diffusion matrix A for weak formulation.
        a (ngsolve.comp.BilinearForm): Gridfunction of order specified by attribute "order" and
                                           lifes in FEspace. Discrete function representing the previous
                                           time step during computation.
                                           \rho_0 = uvold.components[0], c =  uvold.components[1].

        """
        logger.info("Defining weak formulation")
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

    ############################################################################
    ########################## Run and solve system  ###########################
    ############################################################################
    def run_experiment(self):
        """ Run experiment using all class variables defined for the experiment.
            Solution of last time computed successfully is stored in class variable gfucalc.
            successfully computed timesteps are also saved using the ngsolve export for gridfunction.
            If there already saved files, tries to continue the computation from last time step onwards.

            Starts netgen and draws rho if visualoutput_solver is True.

            Creates data folder for experiment and the following files:
            time_list.csv,
            log file with current timestamp.

        """
        logger.info("Starting run_experiment method.")
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
        expri_data = self.expri_data

        if visualoutput_solver == True:
            import netgen.gui

        simulations_folder = self.simulations_folder
        path = self.path
        filedir = self.filedir

        print('simulations saved at: '+path+experiment_full_name)
        try:
            os.mkdir(path+experiment_full_name)
        except Exception as exception:
            logger.error('Simulation folder already existed', exc_info=True)
        ############################################################################
        ############### Continuation of old experiments  ###########################
        ############################################################################
        # checks if experiment data folder already exists, and some steps are already performed.
        cont_switch = folder_checker(experiment_full_name,path)
        logger.debug('cont_switch variable: '+str(cont_switch))
        if cont_switch == True:
            uvold = self.uvold
            old_timestepping_last_time, uvold = continuation_time(experiment_full_name,path,FEspace,dt,uvold)
            if old_timestepping_last_time > 0:
                startingtimestep = int(old_timestepping_last_time/dt)+1
                print(startingtimestep)
            else:
                startingtimestep=0
        elif cont_switch == False:
            uvold = self.uvold
            uvold.Save(filedir+'{}/{}_time_{}'.format(experiment_full_name,experiment_full_name,0))
            startingtimestep = 0
        else:
            print('Something went wrong in the continuation time attempt.')
            logger.debug('something went wrong in the continuation time attempt.', exc_info=True)

        logger.info('Starting time step: '+str(startingtimestep))
        ############################################################################
        ####################### Solve and plot in netgen  ##########################
        ############################################################################

        # computes mass of initial data.
        mass = (Integrate(uvold.components[0],mesh),Integrate(uvold.components[1],mesh))
        gfuL2.Set(uvold.components[0])
        # computes min and max as bad approximation of L^\infty norm.
        linfmin_l2 = min(gfuL2.vec)
        linfmax_l2 = max(gfuL2.vec)
        # Number of time steps.
        nsteps = int(floor(T/dt))
        logger.info('Initial mass: '+str(mass))
        logger.info('Number of time steps: ' + str(nsteps))
        print('Number of time steps: ', nsteps)
        print('Total mass of initial data: ', mass)
        # Plots rho and c, use interface to switch between them.
        if visualoutput_solver == True:
            logger.debug('Plotting initial data')
            gfucalc.vec.data = uvold.vec
            Draw(gfucalc.components[1],mesh, 'c')
            Draw(gfucalc.components[0],mesh, name = 'rho')
            Redraw(True)
        else:
            gfucalc.vec.data = uvold.vec

        # Creates time_list.csv if not continuing computation.
        logger.debug('Removing last line in time_list.csv')
        if cont_switch == False:
            time_list_dict = {'dt': dt, 'timestep': 0,'filename':simulations_folder+experiment_full_name+'/'+experiment_full_name+'_time_0',\
                              'linfymin': round(linfmin_l2,3),'linfmax':  round(linfmax_l2,3),'mass1': round(mass[0],3),'mass2':round(mass[1],3)}

            with open(path+'{}/time_list.csv'.format(experiment_full_name),'w') as csvfile:
                # fieldnames = ['dt','timestep','filename','linfymin','linfmax','mass1','mass2']
                # listwriter = csv.DictWriter(csvfile, fieldnames = fieldnames)
                # listwriter.writerow(time_list_dict)
                pass
        else:
            # deletes last time step of time_list as it was probably not valid.
            with open(path+'{}/time_list.csv'.format(experiment_full_name),'r+') as csvfile:
                fieldnames = ['dt','timestep','filename','linfymin','linfmax','mass1','mass2']
                lines = csvfile.readlines()
                csvfile.truncate(0)
                lines = lines[:-1]
                if len(lines) > 1:
                    lines[-1] = lines[-1].strip()
                csvfile.seek(0)
                for line in lines:
                    csvfile.write(line)

        # creates log file using log_file_creater function of solver.py.
        logger.debug('creates log file using log_file_creater function of solver.py.')
        log_file_name = log_file_creator(path,expri_data)
        sleep(1)
        # creating parameter list for experiment.
        paramlisto = paramlist(path,experiment_full_name,mass,alpha,epsilon, \
                               delta,dt,meshsize,order,T,geometry,ini_data_str)
        # Sets number of threads for parallel computing.
        # Use how many cores/threads you have, as always.
        SetNumThreads(4)
        sleep(1)
        # Actually run the solver of solver.py.
        try:
            logger.info('Start run method from solver.')
            self.gfucalc,endtime = run(path,experiment_full_name,uvold,gfucalc,gfuL2,weakform, \
                                       mesh,dt,nsteps,paramlisto,startingtimestep,visualoutput_solver,log_file_name)
            print('Experiment: {}, done at {}.'.format(experiment_full_name,endtime))
            logger.info('Experiment: {}, done at {}.'.format(experiment_full_name,endtime))
        except Exception as e:
            print(str(e))
            with open(log_file_name,'a') as f:
                f.write(str(e)+';\n')
            logger.info('Some error occured during the run method of solver.py.', exc_info=True)
        finally:
            pass

    ############################################################################
    ########################## Plot functions ##################################
    ############################################################################

    def plotseries_and_more(self, netgenrender = True, vtkexport = False, make_linftyplot = False, \
                            scale_factor = 1, plot_modulo = 1,deform_param = 1):
        """ Plots simulation data in netgen, generates L^\infty plot over time and
            generates vtk output.

        Parameter
        ----------

        experiment_name (str): used in all export or saving file.
        alpha, delta (float): Free parameters to
        epsilon (float): Parameter to switch from parabolic-parabolic to parabolic-elliptic modelself.
        dt (float): time step size.
        T (float): end time of the simulation.
        meshsize (float): mesh size.
        order (float): order of the functions in the finite element space.
        geometry (str): geometry specification as string.
        visualoutput_solver (boolean): True starts netgen for plotting during computations.


        """
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
        path_to_expri = path+experiment_full_name+'/'
        try:
            os.mkdir(path_to_expri+'vtk')
        except Exception as FileExistsError:
            try:
                shutil.rmtree(path_to_expri+'vtk')
            except OSError as e:
                print ("Error: %s - %s." % (e.filename, e.strerror))

            os.mkdir(path_to_expri+'vtk')

        timelist = reorder_timesteps(experiment_full_name,path_to_expri)
        paramdict = opencsvfile(experiment_full_name,path_to_expri)
        Linftymin = []
        Linftymax = []
        t = []
        k = 0

        for i in timelist:
            if i[0] == 0.0:
                gfucalc.Load(str(i[1]))
                gfuL2.Set(gfucalc.components[0])
                linfmin_l2 = min(gfuL2.vec)
                linfmax_l2 = max(gfuL2.vec)
                Linftymin.append(linfmin_l2)
                Linftymax.append(linfmax_l2)

                if netgenrender == True:
                    Draw(gfucalc.components[0],mesh,'density')
                    # visoptions.scalfunction="density"
                    # visoptions.vecfunction = "None"
                    # visoptions.scaledeform1 = scale_factor
                    # visoptions.deformation = deform_param
                    Redraw(True)

                inputparam = input('Press "y" to start the simulations, and "n" to return to python console > ')
                if inputparam == 'n':
                    break

                t.append(i[0])
                if vtkexport == True:
                    vtk = VTKOutput(ma=mesh,coefs=[gfucalc.components[0],gfucalc.components[1]],names=['rho','c'],filename=path_to_expri+'vtk/{}_{}'.format(experiment_full_name,k),subdivision=2)
                    vtk.Do()
            else:
                if k % plot_modulo == 0:
                    print(Integrate(gfucalc.components[0],mesh))
                    gfucalc.Load(str(i[1]))
                    gfuL2.Set(gfucalc.components[0])
                    linfmin_l2 = min(gfuL2.vec)
                    linfmax_l2 = max(gfuL2.vec)
                    Linftymin.append(linfmin_l2)
                    Linftymax.append(linfmax_l2)
                    t.append(float(i[0]))
                    if netgenrender == True:
                        Redraw(True)
                        Draw(gfucalc.components[0],mesh,'density')
                        # visoptions.scalfunction="density"


                    if vtkexport == True:
                        vtk = VTKOutput(ma=mesh,coefs=[gfucalc.components[0],gfucalc.components[1]],names=['rho','c'],filename=path_to_expri+'vtk/{}_{}'.format(experiment_full_name,k),subdivision=2)
                        vtk.Do()
                        print('Timestep: {}'.format(i))

                    print(i[0])

            k += 1


        print('Timesteps: {}'.format(t))
        print('L^infty values: {}'.format(Linftymax))
        fig = plt.figure(1)
        plt.subplot(211)
        plt.plot(t,Linftymax)
        plt.title("Maximumvalues vs Time")
        plt.xlabel("t")
        plt.ylabel("Maximumvalue of rho")
        plt.subplot(212)
        plt.plot(t,Linftymin)
        plt.title("Minimumvalue vs Time")
        plt.xlabel("t")
        plt.ylabel("Minimumvalue of rho")
        z = np.array([t,Linftymax])
        z= np.transpose(z)
        lasttime_plotted = round(T*dt,5)
        if make_linftyplot == True:
            try:
                os.mkdir(path+'Experiments_Linftyplots')
            except Exception as FileExistsError:
                try:
                    shutil.rmtree(path+'Experiments_Linftyplots')
                except OSError as e:
                    print ("Error: %s - %s." % (e.filename, e.strerror))

                os.mkdir(path+'Experiments_Linftyplots')
            print('Exporting L^infty plots to png and dat.')
            df = pd.DataFrame(z, columns = ['t','f(t)'])
            df.to_csv(path+'Experiments_Linftyplots/{}_alpha_{}_delta_{}.dat'.format(experiment_full_name,alpha,delta), index=None, header=None)
            plt.savefig(path+'Experiments_Linftyplots/{}_endtime_{}.png'.format(experiment_full_name,i[0]))
        plt.close(fig)

    def plot_and_export(self, at_times = [0], sol_component = 0,netgenrender = True, \
                        matplotlib_export = False, vtkexport = False, scale_factor = 1, \
                        plot_modulo = 1,deform_param = 1):
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
        path_to_expri = path+experiment_full_name+'/'

        if netgenrender == True:
            import netgen.gui

        print(path_to_expri)
        print(experiment_full_name)

        for i in at_times:
            print(i)
            print(path_to_expri+'{}_time_{}'.format(experiment_full_name,i))
            gfucalc.Load(path_to_expri+'{}_time_{}'.format(experiment_full_name,i))
            if netgenrender == True:
                if sol_component == 0:
                    Draw(gfucalc.components[1],mesh,'u2')
                    Draw(gfucalc.components[0],mesh,'u1')
                else:
                    Draw(gfucalc.components[0],mesh,'u1')
                    Draw(gfucalc.components[1],mesh,'u2')

                # visoptions.scalfunction="density"
                visoptions.vecfunction = "None"
                visoptions.scaledeform1 = 0.001
                visoptions.deformation = 1

            gfuL2.Set(gfucalc.components[sol_component])
            linfmin = min(gfuL2.vec)
            linfmax = max(gfuL2.vec)
            x = []
            y = []
            z = []
            trigs = []
            gfuL2.Set(gfucalc.components[sol_component])
            linfmin_l2 = min(gfuL2.vec)
            linfmax_l2 = max(gfuL2.vec)

            if matplotlib_export == True:
                try:
                    os.mkdir(path_to_expri+'plots')
                except Exception as FileExistsError:
                    pass
                # store vertex coorinates in x/y and function value in z array
                for v in mesh.vertices:
                    p = v.point
                    x.append(p[0])
                    y.append(p[1])
                    z.append(gfucalc.components[sol_component].vec[v.nr])


                # triangulation
                for e in mesh.Elements():
                    trigs.append( [v.nr for v in e.vertices] )

                fig = plt.figure()
                ax = fig.gca(projection='3d')
                ax.grid(True)
                fig.patch.set_facecolor('white')
                surf = ax.plot_trisurf(x, y, z, triangles=trigs, antialiased=False, cmap=cm.coolwarm, linewidth=0, shade=False)
                # Activate for rectangle experiment!!
                # ax.set_yticks([0, 0.1,0.2,0.3,0.4])
                # plt.ylim(0, 0.4)
                plt.xlabel('x')
                plt.ylabel('y')
                #
                # fig.colorbar(surf, shrink=0.5, aspect=10)
                # surf.set_edgecolor(surf.to_rgba(surf._A)) # probably transparent edges
                # surf.set_edgecolor('black')
                # x = surf.get_edgecolors()
                # print(x)

                plt.savefig(path_to_expri+'plots'+'/{}_time_{}.eps'.format(experiment_full_name,round(float(i),5)))

                print('Minimum: {}'.format(linfmin_l2))
                print('Maximum: {}'.format(linfmax_l2))
                plt.draw()
                plt.waitforbuttonpress(0)
                plt.close(fig)
