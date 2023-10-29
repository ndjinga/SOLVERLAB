#!/usr/bin/env python
# -*- coding: utf-8 -*-

# %% LICENSE_SALOME_CEA_BEGIN
# see PACKAGESPY/LICENCE file
# %% LICENSE_END


import os
import glob
import unittest
import xyzpy.test.modelSimpleTest as MXYZ
import xyzpy.controllerXyz as CXYZ
import xml.etree.ElementTree as ET
from xyzpy.guiXyz.treeXmlXyz import TreeXmlXyz, TreeXmlXyzItem
import xyzpy.utilsXyz as UXYZ

from salomepy.onceQApplication import OnceQApplication
import xyzpy.intFloatListXyz as IFLX #append factory classes
import xyzpy.utilsXyz as UXYZ
import xyzpy.loggingXyz as LOG
logger = LOG.getLogger()

setTimer = True
deltaTime = 500
withShow = True
verbose = False
PUSHLEVEL = "CRITICAL"

class TestCase(unittest.TestCase):
  def launchTimer(self, app, wid):
    timer = QtCore.QTimer();
    timer.timeout.connect(wid.close)
    if setTimer: timer.start(deltaTime)
    app.exec_()

  def test_010(self):
    a = CXYZ.ControllerXyz()
    self.assertEqual("ControllerXyz[" in  a.objectName(), True)
    LOG.pushLevel(PUSHLEVEL)
    self.assertEqual(a.getModel(), None)
    LOG.popLevel()
    self.assertEqual(a.getRequest("hello...").typeRequest, "hello...")
    
  def test_015(self):
    # pattern MCV with no View with services signals from API controller
    C = CXYZ.ControllerXyz()
    M = MXYZ.ModelSimpleTest()
    C.setModel(M)
    self.assertEqual(C._model, M)
    self.assertEqual(M.getController(), C)
    self.assertEqual(M.lastReceiveRequest(), None)
    
  def test_020(self):
    # pattern MVC with every services signals from API controller
    app = OnceQApplication()
    C = CXYZ.ControllerXyz()
    M = MXYZ.ModelSimpleTest()
    C.setModel(M)
    V1 = TreeXmlXyz()
    C.setView(V1)
    
    self.assertNotEqual(C.getModel(), M) #getModel send xml copy of model
    self.assertEqual(M.getController(), C)
    self.assertEqual(V1.getController(), C)
    self.assertEqual(C.getViews(), [V1])

    self.assertRaises(Exception , C.setView, V1)
    V2 = TreeXmlXyz()
    C.setView(V2)
    self.assertEqual(V2.getController(), C)
    self.assertEqual(C.getViews(), [V1, V2]) #alphabetical sort on Views objectName
    
    self.assertNotEqual(C.getModel(), M)
    self.assertEqual(M.lastReceiveRequest(), None)
    self.assertEqual(V1.lastReceiveRequest(), None)
    self.assertEqual(V2.lastReceiveRequest(), None)

    aRequest = C.getRequest("hello from V1")
    self.assertEqual(V1.getController().sendRequestToController(aRequest), True)
    self.assertNotEqual(M.lastReceiveRequest(), None)
    self.assertEqual("hello from V1</typeRequest" in M.lastReceiveRequest(), True)
    
    aRequest = C.getRequest("hello from V2")
    self.assertEqual(V2.getController().sendRequestToController(aRequest), True)
    self.assertNotEqual(M.lastReceiveRequest(), None)
    self.assertEqual("hello from V2</typeRequest" in M.lastReceiveRequest(), True)
    
    aRequest = C.getRequest("hello from M")
    self.assertEqual(M.getController().sendRequestToViews(aRequest), True)
    self.assertNotEqual(M.lastReceiveRequest(), None)
    self.assertEqual("hello from V2</typeRequest" in M.lastReceiveRequest(), True)
    self.assertEqual("hello from M</typeRequest" in V1.lastReceiveRequest(), True)
    self.assertEqual("hello from M</typeRequest" in V2.lastReceiveRequest(), True)
    
    if verbose: print("\n********test_020_1\n")
    
    #test send to model, and reply to Views
    aRequest = V1.getController().getRequest("TestRequestFromView hello from V1") #test V->C->M->C->V
    self.assertNotEqual(M.lastReceiveRequest(), None)
    V1.getController().sendRequestToController(aRequest)
    self.assertNotEqual(M.lastReceiveRequest(), None)
    self.assertEqual("TestRequestFromView hello from V1</typeRequest" in M.lastReceiveRequest(), True)
    self.assertEqual("TestRequestFromView hello from V1</typeRequest" in V1.lastReceiveRequest(), True)
    self.assertEqual("TestRequestFromView hello from V1</typeRequest" in V2.lastReceiveRequest(), True)
    aRequest = V2.getController().getRequest("TestRequestFromView helloMore from V2") #test V->C->M->C->V
    V2.getController().sendRequestToController(aRequest)
    self.assertNotEqual(M.lastReceiveRequest(), None)
    self.assertEqual("TestRequestFromView helloMore from V2</typeRequest" in M.lastReceiveRequest(), True)
    self.assertEqual("TestRequestFromView helloMore from V2</typeRequest" in V1.lastReceiveRequest(), True)
    self.assertEqual("TestRequestFromView helloMore from V2</typeRequest" in V2.lastReceiveRequest(), True)
    
 
  def test_030(self):
    # pattern MVC with every services signals from API controller
    app = OnceQApplication()
    C = CXYZ.ControllerXyz()
    #Request
    aRequest = C.getRequest("hello...")
    #theorically done by aView: aView.getController().sendRequestToController(aRequest)
    LOG.pushLevel(LOG.levels["CRITICAL"]+1)
    self.assertEqual(C.sendRequestToController(aRequest), False) #no model, but never mind, warning
    LOG.popLevel()
    M = MXYZ.ModelSimpleTest()
    C.setModel(M)
    self.assertEqual(C.sendRequestToController(aRequest), True) #a model, ok


if __name__ == '__main__':
  setTimer = False #True #False #user call and wait
  from salomepy import qMainWindowForLog
  qMainWindowForLog.verboseEvent = not setTimer
  from xyzpy.guiXyz import treeXmlXyz
  #treeXmlXyz.verboseEvent = not setTimer
  unittest.main()
  pass
