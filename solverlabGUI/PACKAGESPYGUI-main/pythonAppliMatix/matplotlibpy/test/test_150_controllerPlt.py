#!/usr/bin/env python
# -*- coding: utf-8 -*-

# %% LICENSE_SALOME_CEA_BEGIN
# see PACKAGESPY/LICENCE file
# %% LICENSE_END


import os
import glob
import unittest
from PyQt5 import QtGui, QtCore, QtWidgets

from matplotlibpy.controllerPlt import ControllerPlt
import xyzpy.utilsXyz as UXYZ

from salomepy.onceQApplication import OnceQApplication
import xyzpy.intFloatListXyz as IFLX #append factory classes
import xyzpy.loggingXyz as LOG
logger = LOG.getLogger()


setTimer = True
deltaTime = 500
withShow = True
verbose = False
PUSHLEVEL = "CRITICAL"
testDir = os.path.realpath(
          os.path.split(
          os.path.realpath(__file__))[0])

class TestCase(unittest.TestCase):
  def launchTimer(self, app, wid):
    timer = QtCore.QTimer();
    timer.timeout.connect(wid.close)
    if setTimer: timer.start(deltaTime)
    app.exec_()

  def test_001(self):
    LOG.pushLevel(PUSHLEVEL)
    app = OnceQApplication()
    desktop = QtWidgets.QMainWindow()
    desktop.resize(1000, 600)
    ctrl = ControllerPlt(desktop = desktop)
    desktop.setWindowTitle(ctrl.objectName())
    #aModel = UXYZ.fromFileXml(os.path.join(testDir,"aMatplotlibPltExample.xml"))
    #print "aModel", aModel
    #ctrl.setModel(aModel)
    desktop.show()
    self.launchTimer(app, desktop)
    LOG.popLevel()

  def test_010(self):
    app = OnceQApplication()
    LOG.pushLevel(PUSHLEVEL)
    a = ControllerPlt()
    self.assertEqual("ControllerPlt[" in  a.objectName(), True)
    self.assertEqual(a.getModel(), None)
    self.assertEqual(a.getRequest("hello...").typeRequest, "hello...")
    LOG.popLevel()

  def test_015(self):
    app = OnceQApplication()
    a = ControllerPlt()
    fig1 = 1 #no good type but it do not care...
    a.addNewFigure(fig1) #set default name for a plot
    fig2 = 2
    a.addNewFigure(fig2)
    a.addNewFigure(fig2)
    res = [i.name for i in a._model]
    self.assertEqual(len(res), 3)
    a.deletePlt(res[1])
    res2 = [i.name for i in a._model]
    self.assertEqual(len(res2), 2)
    self.assertEqual(res2, [res[0], res[2]])

  def test_100(self):
    LOG.pushLevel(PUSHLEVEL)
    app1 = ControllerPlt()
    self.assertNotEqual(app1, None)
    app2 = ControllerPlt()
    self.assertNotEqual(app1, app2)
    self.assertNotEqual(id(app1), id(app2))
    LOG.popLevel()

  def test_110(self):
    LOG.pushLevel(PUSHLEVEL)
    ctrl = ControllerPlt()
    aRequest = ctrl.getRequest("LaunchAllTestsAction")
    ctrl.sendRequestToController(aRequest, verbose=False) # no model
    LOG.popLevel()


if __name__ == '__main__':
  verbose = False
  if verbose: PUSHLEVEL = "INFO"
  setTimer = False #False #user call and wait
  unittest.main()
  pass
