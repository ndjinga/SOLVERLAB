#!/usr/bin/env python
# -*- coding: utf-8 -*-

# %% LICENSE_SALOME_CEA_BEGIN
# see PACKAGESPY/LICENCE file
# %% LICENSE_END

"""\
module managing paraview display on miscellanous types of files
see ./example/paraviewDisplayExample.py
"""

import os
from PyQt5 import QtGui,QtCore, QtWidgets
import salomepy.xsalomesession as XSS

import pvsimple as PV


def paraviewDisplayMicrostructureVtk(self):
    """display from self is nameFile.vtk"""
    
    try:
      desktop = self.getController().getDesktop()
    except:
      #desktop = None
      desktop = XSS.getDesktop()
      
    PV.ShowParaviewView()
    fileName = str(self)
    _, name = os.path.split(fileName)
    
    if not os.path.isfile(fileName):
      QtWidgets.QMessageBox.warning(desktop, "warning", 
          "Display microstructure vtk data: not a file \n'%s'" % fileName)
      return False

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

    PV.Render()
    XSS.updateObjBrowser()
    return True

  
