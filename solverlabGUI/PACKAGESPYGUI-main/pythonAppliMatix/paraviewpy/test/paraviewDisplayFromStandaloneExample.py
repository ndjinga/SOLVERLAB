#!/usr/bin/env python
# -*- coding: utf-8 -*-

# %% LICENSE_SALOME_CEA_BEGIN
# see PACKAGESPY/LICENCE file
# %% LICENSE_END

"""\
test module managing paraview display on miscellanous types of files
"""

import os
import pprint as PP

"""
#export OMNIORB_CONFIG='/tmp/.omniORB_wambeke_is210329_2813.cfg'
OMNIORB_CONFIG = os.getenv("OMNIORB_CONFIG")
if OMNIORB_CONFIG == None:
  OMNIORB_CONFIG = '/tmp/.omniORB_wambeke_is210329_2814.cfg'
  os.environ["OMNIORB_CONFIG"] = OMNIORB_CONFIG

print "OMNIORB_CONFIG='%s'" % os.getenv("OMNIORB_CONFIG")

#import paraviewpy.pvsimple_cv as PV
"""

import os
import time
import salomepy.xsalomesession as XSS
XSS.set_OMNIORB_CONFIG()

import pvsimple as PV

def essai():
    testDir = os.path.join(os.path.split(os.path.realpath(__file__))[0])
    aVtkFile = os.path.join(testDir, "voronoi_10grains_voxelized.vtk")
    
    #PV.ShowParaviewView()
    fileName = str(aVtkFile)
    _, name = os.path.split(fileName)
    
    #### disable automatic camera reset on 'Show'
    PV._DisableFirstRenderCameraReset()
    # create a new 'Legacy VTK Reader'
    zone_vtk = PV.LegacyVTKReader(FileNames=[fileName])
    PV.RenameSource(name)
    # get animation scene
    animationScene1 = PV.GetAnimationScene()
    # update animation scene based on data timesteps
    animationScene1.UpdateAnimationUsingDataTimeSteps()
    # get active view
    renderView1 = PV.GetActiveViewOrCreate('RenderView')
    # uncomment following to set a specific view size
    # renderView1.ViewSize = [697, 713]
    # show data in view
    zone_vtkDisplay = PV.Show(zone_vtk, renderView1)
    # reset view to fit data
    #renderView1.ResetCamera()
    # create a new 'Threshold'
    threshold1 = PV.Threshold(Input=zone_vtk)
    # show data in view
    threshold1Display = PV.Show(threshold1, renderView1)
    # show color bar/color legend
    threshold1Display.SetScalarBarVisibility(renderView1, True)
    # get color transfer function/color map for 'MaterialId'
    materialIdLUT = PV.GetColorTransferFunction('MaterialId')
    # get opacity transfer function/opacity map for 'MaterialId'
    materialIdPWF = PV.GetOpacityTransferFunction('MaterialId')
    # rescale color and/or opacity maps used to exactly fit the current data range
    threshold1Display.RescaleTransferFunctionToDataRange(False)
    
    camera = PV.GetActiveCamera()
    #pos = camera.GetPosition()
    #print "pos",pos
    camera.SetPosition(-1000, 1000, -2000) #point of view far away
    renderView1.ResetCamera() #fit zoom keeping point of view
    
    if PV.__getFromGUI() != 0: #TODO 161026 solve with adrien PV.Render() create window
      PV.Render()
      XSS.updateObjBrowser()
    
    return True

print("mode __name__", __name__)
XSS.set_OMNIORB_CONFIG()
OMNIORB_CONFIG = os.getenv("OMNIORB_CONFIG")
print("OMNIORB_CONFIG", OMNIORB_CONFIG)
print("oops")
  
if __name__ == '__main__':
  """
  Identify if we are running inside SALOME's embedded interpreter.
  @return a value strictly greater than 0 if we are in SALOME's embedded interpreter
  @return 2 if we are in Salome embedded Python console.
  """
  print("PV.__getFromGUI()=%s : TODO 161026 solve with adrien PV.Render() create window" % PV.__getFromGUI())
  essai()
  time.sleep(20) #wait a time to see PV.Render() create window if done
  """
  print "dir():\n%s" % PP.pformat(dir())
  print "dir(PV):\n%s" % PP.pformat(dir(PV))
  print PP.pformat(dir(PV.GetRenderView()))
  print PV.GetRenderViews()
  print PV.GetRenderView()
  print PV.GetRenderView().ViewSize
  """

    
