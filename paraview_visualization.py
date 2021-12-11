#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

N=40

out_dir='/scratch/ws/1/haja565a-haja565a/master_thesis/output'



# create a new 'PVD Reader'
file0_pvd = PVDReader(FileName=out_dir + '/phase'+ '_p0_.pvd')


# get active view
renderView = GetActiveViewOrCreate('RenderView')

# show data in view
fileDisplay = Show(file0_pvd, renderView)
fileDisplay.SetScaleArray = ['POINTS', 'phi']

# set scalar coloring
ColorBy(fileDisplay, ('POINTS', 'phi'))

# get animation scene
animationScene = GetAnimationScene()

# update animation scene based on data timesteps
animationScene.UpdateAnimationUsingDataTimeSteps()

for i in range(N):
  #if i != 18 and i != 158:
    #continue
  
  # create a new 'PVD Reader'
  # file_pvd = PVDReader(FileName=out_dir + '/phase' + postfix + '_p' + str(i) + '_.pvd')  
  file_pvd = PVDReader(FileName=out_dir + '/phase' +  '_p' + str(i) + '_.pvd')
  # create a new 'Threshold'

  threshold = Threshold(Input=file_pvd, Scalars=['POINTS', 'phi'], ThresholdRange=[-0.95, 1.2])

  # show data in view
  thresholdDisplay = Show(threshold, renderView)
  thresholdDisplay.SetScaleArray = ['POINTS', 'phi']

  # set scalar coloring
  ColorBy(thresholdDisplay, ('POINTS', 'phi'))
  
 # print ((i*100.0)/N,'%')	


# reset view to fit data
renderView.ResetCamera()

#changing interaction mode based on data extents
renderView.InteractionMode = '2D'
renderView.CameraPosition = [50.0, 50.0, 10000.0]
renderView.CameraFocalPoint = [50.0, 50.0, 0.0]
