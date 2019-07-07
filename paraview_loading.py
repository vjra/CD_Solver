from paraview.simple import *
import csv
import os


experiment_name ='Hitt_Expri_42_alpha_1_epsilon_0_delta_1'
# screenshot_foldername = 'end_time_screenshots'
# Scale_Factor = 0.5
# timestep = -1
# # makes video and ignores timestep variable
# make_video = False
# interactive_mode = False

path = '/mnt/data/simulations/paper/'

def log_file_creator(foldername,experiment_list):
    today = datetime.now()
    today_string = '{}_{}_{}_{}'.format(today.year,today.month,today.day,today.hour)
    logfilename = path+foldername+'/{}_errorlog.csv'.format(today_string)
    try:
        with open(logfilename,"w") as f:
            f.write("##############LOG############\n")
            for expridata in experiment_list:
                f.write(str(expridata)+';\n')
    except Exception as FileExistsError:
        pass

    return logfilename

def opencsvfile(name):

    with open(path+'/simulations/none/{}/parameter.csv'.format(name), 'r') as csvfile:
        csv_reader = csv.DictReader(csvfile)
        paramdict = list(csv_reader)

    return paramdict[0]

def load_vtk_data(experiment_name):
    paramdict = opencsvfile(experiment_name)
    T = float(paramdict['T'])
    dt = float(paramdict['dt'])
    maxtimeT = int(T/dt)
    print(T,dt,maxtimeT)
    vtklist = []
    number_of_vtkfiles = 0
    for i in range(0,maxtimeT,1):
        if True == os.path.isfile(path+'/simulations/none/{}/vtk/{}_{}.vtk'.format(experiment_name,experiment_name,i)):
            vtklist.append(path+'simulations/{}/vtk/{}_{}.vtk'.format(experiment_name, experiment_name,i))
            number_of_vtkfiles = i

    print('Number of vtk-files in folder: {}'.format(number_of_vtkfiles))
    return vtklist,T,dt,number_of_vtkfiles




def export_screenshot(experiment_name,Scale_Factor, timestep, make_video,foldername,interactive_mode):
    #### disable automatic camera reset on 'Show'
    #paraview.simple._DisableFirstRenderCameraReset()

    ############################ Load saved state from PV ############################

    # servermanager.LoadState('test_state.pvsm')

    ############################ Load .vtk data ############################

    vtklist,T,dt,number_of_vtkfiles = load_vtk_data(experiment_name)
    screenshot_time = float(number_of_vtkfiles*dt)


    if make_video == True:
        reader = OpenDataFile(vtklist)
    else:
        reader = OpenDataFile(vtklist[timestep])


    ############################ make nice rendering with paraview ############################

    # set activ source
    SetActiveSource(reader)
    # get active view
    renderView1 = GetActiveViewOrCreate('RenderView')
    pd = reader.PointData

    # for n in range(pd.GetNumberOfArrays()):
    #     print "Range for array ", pd.GetArray(n).GetName(), " ", pd.GetArray(n).GetRange()
    rho_data_range = pd.GetArray(1).GetRange()
    print("Maxvalue: {}, Minvalue: {}".format(int(rho_data_range[1]),int(rho_data_range[0])))
    negativevalue_appears = False
    if float(rho_data_range[0]) <0:
        negativevalue_appears = True
    if Scale_Factor == True:
        maxrho_scale = 2/float(int(rho_data_range[1]))
    else:
        maxrho_scale = 1

    # update pipeline
    UpdatePipeline()
    ksdisplay1 = Show(reader,renderView1)
    #ColorBy(ksdisplay1, ('POINTS', 'rho'))
    #ksdisplay1.RescaleTransferFunctionToDataRange(True)
    #ksdisplay1.SetScalarBarVisibility(renderView1, True)

    renderView1.ResetCamera()

    # make new pipeline object
    pipe00 = GetActiveSource()
    Hide(reader)

    # make deformation with warp scalar, scaled and with scalar rho
    warping = WarpByScalar(Input = pipe00, Scalars = 'rho', ScaleFactor = maxrho_scale)
    #print(warping.ListProperties())


    warpingdisplay1 = Show(warping,renderView1)
    # color by value of rho
    ColorBy(warpingdisplay1, ('POINTS', 'rho'))
    UpdatePipeline()
    # styling, background white, activate axis grid and label it black,
    #show grid, and scale it such that it fits the scaling of wrapbyscalar,
    renderView1.Background = [1,1,1]
    renderView1.AxesGrid
    renderView1.AxesGrid.XTitle = 'x'
    renderView1.AxesGrid.YTitle = 'y'
    renderView1.AxesGrid.ZTitle = 'rho'
    renderView1.AxesGrid.XTitleColor = [0.0, 0.0, 0.0]
    renderView1.AxesGrid.YTitleColor = [0.0, 0.0, 0.0]
    renderView1.AxesGrid.ZTitleColor = [0.0, 0.0, 0.0]
    renderView1.AxesGrid.Visibility = 1
    renderView1.AxesGrid.ShowGrid = 0
    renderView1.AxesGrid.GridColor = [0.0, 0.0, 0.0]
    renderView1.AxesGrid.DataScale = [1, 1, maxrho_scale]
    renderView1.OrientationAxesVisibility = 0
    # make fixed window size
    renderView1.ViewSize = [1024, 863]
    warpingdisplay1.RescaleTransferFunctionToDataRange(True)
    warpingdisplay1.SetScalarBarVisibility(renderView1, False)
    # creating good camera position and angle
    renderView1.CameraPosition = [1.5, -6, 4]
    renderView1.CameraFocalPoint = [0.0, 0.0, 1.5]
    renderView1.CameraViewUp = [-0.1, 0.5, 1]
    renderView1.CameraParallelScale = 1.2
    # resets the camera, such that the whole plot is visible
    renderView1.ResetCamera()
    # render and show image
    Render()
    # activates interactive mode
    if interactive_mode == True:
        Interact()
    # after interaction, we do again a reset of the Camera
    # such that all of the plot is visible before we make a screenshot
    renderView1.ResetCamera()


    ############################ Play animation ############################
    # animationScene1 = GetAnimationScene()
    # animationScene1.Play()
    ############################ Save animation ############################
    #save animation images/movie
    if make_video == True:
        WriteAnimation(path+'/videos/{}'.format(experiment_name)+'.avi', Magnification=1, FrameRate=25.0, Compression=True)


    # export scene for pdf, eps and svg
    # ExportView('/home/oleingan/Nextcloud/KS2/ngsolve/repos/KS_ngsolve/nr2.eps', view=renderView1, Plottitle='ParaView GL2PS Export',
    #     Compressoutputfile=0,
    #     Drawbackground=0,
    #     Cullhiddenprimitives=1,
    #     Linewidthscalingfactor=0.714,
    #     Pointsizescalingfactor=0.714,
    #     GL2PSdepthsortmethod='Simple sorting (fast, good)',
    #     #GL2PSdepthsortmethod='BSP sorting (slow, best)',
    #     Rasterize3Dgeometry=0,
    #     Dontrasterizecubeaxes=0,
    #     Rendertextaspaths=0)
    ############################ Save screenshot/scene ############################
    # # save screenshot for png, jpeg, etc.
    if negativevalue_appears == True:
        SaveScreenshot(path+'{}/NEG_{}_at_time_{}.png'.format(foldername,experiment_name,screenshot_time), quality=100, view=renderView1)
    else:
        SaveScreenshot(path+'{}/{}_at_time_{}.png'.format(foldername,experiment_name,screenshot_time), quality=100, view=renderView1)

    Hide()



# export_screenshot(experiment_name,Scale_Factor, timestep, make_video,screenshot_foldername,interactive_mode)
