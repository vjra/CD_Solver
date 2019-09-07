from solver import *
from ngsolve import *
from plots_and_more_lib import *
# import netgen.gui
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

        cont_switch = folder_checker(experiment_full_name,path)
        print('Contswitch variable: '+str(cont_switch))
        if cont_switch == True:
            old_timestepping_last_time, uvold = continuation_time(experiment_full_name,path,FEspace,dt)
            startingtimestep = int(old_timestepping_last_time/dt)+1
        elif cont_switch == False:
            uvold = self.uvold
            uvold.Save(filedir+'{}/{}_time_{}'.format(experiment_full_name,experiment_full_name,0))
            startingtimestep = 0
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
        if cont_switch == False:
            time_list_dict = {'dt': dt, 'timestep': 0,'filename':simulations_folder+experiment_full_name+'/'+experiment_full_name+'_time_0',\
                              'linfymin': round(linfmin_l2,3),'linfmax':  round(linfmax_l2,3),'mass1': round(mass[0],3),'mass2':round(mass[1],3)}
            with open(path+'{}/time_list.csv'.format(experiment_full_name),'w') as csvfile:
                # fieldnames = ['dt','timestep','filename','linfymin','linfmax','mass1','mass2']
                # listwriter = csv.DictWriter(csvfile, fieldnames = fieldnames)
                # listwriter.writerow(time_list_dict)
                pass
        else:
            with open(path+'{}/time_list.csv'.format(experiment_full_name),'r+') as csvfile:
                fieldnames = ['dt','timestep','filename','linfymin','linfmax','mass1','mass2']
                lines = csvfile.readlines()
                lines = lines[:-1]
                lines[-1] = lines[-1].strip()
                csvfile.seek(0)
                print(lines)
                for line in lines:
                    csvfile.write(line)


        sleep(1)
        log_file_name = log_file_creator(simulations_folder,experiment_full_name)
        paramlisto = paramlist(path,experiment_full_name,mass,alpha,epsilon,delta,dt,meshsize,order,T,geometry,ini_data_str)
        SetNumThreads(8)
        sleep(1)
        try:
            self.gfucalc,endtime = run(path,experiment_full_name,uvold,gfucalc,gfuL2,weakform,mesh,dt,nsteps,paramlisto,startingtimestep,visualoutput_solver)
            print('Experiment: {}, done'.format(experiment_full_name))
        except Exception as e:
            print(str(e))
            with open(log_file_name,'a') as f:
                f.write(str(e)+';\n')
        finally:
            pass

    def plotseries_and_more(self, netgenrender = True, vtkoutput = False, make_linftyplot = False,scale_factor = 1, plot_modulo = 1,deform_param = 1):
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

        try:
            os.mkdir(path_to_expri+'vtk')
        except Exception as FileExistsError:
            try:
                shutil.rmtree(path_to_expri+'vtk')
            except OSError as e:
                print ("Error: %s - %s." % (e.filename, e.strerror))

            os.mkdir(path_to_expri+'vtk')

        timelist = reorder_timesteps(experiment_full_name,path_to_expri)
        print(timelist)
        paramdict = opencsvfile(experiment_full_name,path_to_expri)
        Linftymin = []
        Linftymax = []
        t = []
        k = 0

        for i in timelist:
            if i[0] == 0.0:
                gfucalc.Load(str(i[1]))
                gfuL2.Set(gfucalc.components[0])
                # pdb.set_trace()
                x = []
                y = []
                z = []
                trigs = []
                for v in mesh.vertices:
                    p = v.point
                    x.append(p[0])
                    y.append(p[1])
                    z.append(gfucalc.components[0].vec[v.nr])
                linfmax = max(z)
                linfmin = min(z)


                if netgenrender == True:
                    Draw(gfucalc.components[0],mesh,'density')
                    # visoptions.scalfunction="density"
                    visoptions.vecfunction = "None"
                    visoptions.scaledeform1 = scale_factor
                    visoptions.deformation = deform_param

                inputparam = input('Press "y" to start the simulations, and "n" to return to python console > ')
                if inputparam == 'n':
                    break

                Linftymin.append(linfmin)
                Linftymax.append(linfmax)
                t.append(i[0])
                if vtk_export == True:
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


                    # store vertex coorinates in x/y and function value in z array

                    if vtkoutput == True:
                        vtk = VTKOutput(ma=mesh,coefs=[gfucalc.components[0],gfucalc.components[1]],names=['rho','c'],filename=path_to_expri+'vtk/{}_{}'.format(experiment_full_name,k),subdivision=2)
                        vtk.Do()
                        print('Timestep: {}'.format(i))

                    print(i[0])
            k += 1


        print(t)
        print(Linftymax)
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
        # plt.plot(t[0:-2],Linftymax[0:-2])
        # plt.plot(t[0:-2],Linftymin[0:-2])
        #plt.draw()
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

            df = pd.DataFrame(z, columns = ['t','f(t)'])
            df.to_csv(path+'Experiments_Linftyplots/{}_alpha_{}_delta_{}.dat'.format(experiment_full_name,alpha,delta), index=None, header=None)
            plt.savefig(path+'Experiments_Linftyplots/{}_endtime_{}.png'.format(experiment_full_name,i[0]))
        plt.close(fig)

    def plot_and_export(self, at_times = [0], sol_component = 0,netgenrender = True, matplotlib_export = False, vtk_export = False, scale_factor = 1, plot_modulo = 1,deform_param = 1):
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
