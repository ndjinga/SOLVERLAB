#!/usr/bin/env python
# -*- coding: utf-8 -*-

# %% LICENSE_SALOME_CEA_BEGIN
# see PACKAGESPY/LICENCE file
# %% LICENSE_END

"""\
test module managing paraview display on miscellanous types of files
"""

import os
import unittest
import time
import pprint as PP
import threading

try:
  import salomepy.xsalomesession as XSS
  modules = ["GEOM", "SMESH", "MED", "PARAVIS"]
  XSS.getXSalomeSession(showDesktop=True, modules=modules)
except:
  XSS = None
   
try:
  import pvsimple as PV
except:
  PV = None

try:
  import paraviewpy.paraviewDisplay as PVD
except:
  PVD = None
  
setTimer = True
deltaTime = 30
withShow = True
verbose = False
ko = [True] #mutable from TestCase

#print "dir():\n%s" % PP.pformat(dir())

class TestCase(unittest.TestCase):
  
  testDir=os.path.join(os.path.split(os.path.realpath(PVD.__file__))[0], "test")
  aVtkFile = os.path.join(testDir, "voronoi_10grains_voxelized.vtk")

  def launchTimer(self, deltaTime):
    print("time.sleep",deltaTime)
    self.t = threading.Timer(float(deltaTime), self.end)
    self.t.start()
    #time.sleep(deltaTime)
    
  def end(self):
    print("!!!end TestCase!!!")
    
  def test_001(self):
    self.assertNotEqual(XSS, None)
    self.assertNotEqual(PV, None)
    self.assertNotEqual(PVD, None)
    ko[0] = False
    if verbose: print("test_001 ok",not(ko[0]))
    
  def test_100(self):
    if verbose: print("test_100 ok",not(ko[0]))
    if ko[0]: return
    #PVD.paraviewDisplayMicrostructureVtk("/home/wambeke/voronoi_10grains_voxelized.vtk")
    PVD.paraviewDisplayMicrostructureVtk(self.aVtkFile)
    if verbose: print("PVD.paraviewDisplayMicrostructureVtk(%s)" % self.aVtkFile)
    if setTimer: self.launchTimer(deltaTime)


if __name__ == '__main__':
  setTimer = True #True #False #user call and wait
  verbose = True
  print("__main__ of paraviewDisplayTest.py")
  unittest.main()
  pass

 
    
