#!/usr/bin/env python
# -*- coding: utf-8 -*-

# %% LICENSE_SALOME_CEA_BEGIN
# see PACKAGESPY/LICENCE file
# %% LICENSE_END


import os
import sys
import unittest

from PyQt5 import QtCore, QtGui, QtWidgets
import xyzpy.controllerXyz as CXYZ
import xml.etree.ElementTree as ET
from salomepy.onceQApplication import OnceQApplication
from xyzpy.guiXyz.treeXmlXyz import TreeXmlXyz, TreeXmlXyzItem, TreeXmlXyzMainWidget
import xyzpy.intFloatListXyz as IFLX #append factory classes
import xyzpy.utilsXyz as UXYZ
import xyzpy.classFactoryXyz as CLFX
import xyzpy.directoryManagerXyz as DMGX


setTimer = True
deltaTime = 1000
withShow = True
verbose = False
user = os.getenv("USERNAME")

# directory xyzpy
testDir = os.path.realpath(os.path.join(os.path.dirname(__file__), ".."))

class TestCase(unittest.TestCase):

  def root_path():
    return os.path.abspath(os.sep)
  
  def launchTimer(self, app, wid):
    timer = QtCore.QTimer();
    timer.timeout.connect(wid.close)
    if setTimer: timer.start(deltaTime)
    app.exec_()
      
  def test_100(self):
    app = OnceQApplication()
    dirName = os.path.join(testDir)
    self.assertTrue(os.path.isdir(dirName)) #inexisting dir is not a test
    aData = DMGX.DirectoryManagerXyz()
    aData.directoryName = dirName
    aStrData = aData.toStrXml()
    data = ET.fromstring(aStrData)
    treeWid = TreeXmlXyz()
    treeWid.setFromXml(data)
    treeWid.resize(500, 500)
    treeWid.show()
    self.launchTimer(app, treeWid)

  def xtest_200(self):
    dirName = os.path.join(testDir)
    self.assertTrue(os.path.isdir(dirName)) #inexisting dir is not a test
    app = OnceQApplication()
    C = CXYZ.ControllerXyz()
    M = DMGX.DirectoryManagerXyz()
    M.directoryName = dirName
    C.setModel(M)
    self.assertqual(C._model, M)
    self.assertEqual(M.getController(), C)
    self.assertEqual(M.lastReceiveRequest(), None)
    self.launchTimer(app, treeWid)

  def xtest_200(self):
    dirName = os.path.join(testDir)
    self.assertTrue(os.path.isdir(dirName)) #inexisting dir is not a test
    app = OnceQApplication()
    desktop = QtWidgets.QMainWindow()
    desktop.resize(700, 600)
    C = CXYZ.ControllerXyz(desktop = desktop)
    desktop.setWindowTitle("directoryManagerXyz")
    M = DMGX.DirectoryManagerXyz()
    M.directoryName = dirName
    C.setModel(M)
    desktop.show()
    self.launchTimer(app, desktop)

if __name__ == '__main__':
  verbose = True
  setTimer = False #user call and wait
  unittest.main()
  pass
