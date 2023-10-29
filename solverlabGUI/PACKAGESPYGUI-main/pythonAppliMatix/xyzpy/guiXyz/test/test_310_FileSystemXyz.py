#!/usr/bin/env python
# -*- coding: utf-8 -*-

# %% LICENSE_SALOME_CEA_BEGIN
# see PACKAGESPY/LICENCE file
# %% LICENSE_END


"""
>> highlight -h
    USAGE: etc...

>> highlight --list-themes
    Installed themes (located in /usr/share/highlight/themes/):
    acid  etc...

>> highlight --list-langs
    Installed language definitions (located in /usr/share/highlight/langDefs/):
    ABAP/4 etc...
"""


import os
import sys
import unittest

from PyQt5 import QtCore, QtGui, QtWidgets
from salomepy.onceQApplication import OnceQApplication
import subprocess as SP
import mimetypes
import xyzpy.guiXyz.FileSystemXyz as FISY

setTimer = True
deltaTime = 1000
verbose = False
verboseView = False
verboseEvent = False

testDir = os.path.realpath(
          os.path.split(
          os.path.realpath(__file__))[0])



############################################
class TestCase(unittest.TestCase):
  
  def launchTimer(self, app, wid):
    if setTimer: 
      timer = QtCore.QTimer();
      timer.timeout.connect(wid.close)
      timer.start(deltaTime)
    app.exec_()
  
  def xtest_100(self):
    app = OnceQApplication()
    currenDir = os.path.split(testDir)[0] #os.getcwd()
    treeView = QtWidgets.QTreeView()
    fileSystemModel = QtWidgets.QFileSystemModel(treeView)
    fileSystemModel.setReadOnly(True)
    root = fileSystemModel.setRootPath(currenDir)
    treeView.setModel(fileSystemModel)
    treeView.setRootIndex(root)
    treeView.setWindowTitle("test_100")
    treeView.show()
    self.launchTimer(app, treeView)

  def xtest_200(self):
    app = OnceQApplication()
    currenDir = os.path.split(testDir)[0] #os.getcwd()
    treeView = FISY.QfileViewTest()
    fileSystemModel = QtWidgets.QFileSystemModel(treeView)
    fileSystemModel.setReadOnly(True)
    root = fileSystemModel.setRootPath(currenDir)
    treeView.setModel(fileSystemModel)
    treeView.setRootIndex(root)
    treeView.setWindowTitle("test_200")
    treeView.show()
    self.launchTimer(app, treeView)

  def xtest_250(self):
    app = OnceQApplication()
    path = "$PACKAGESPY_ROOT_DIR"
    dialog = FISY.UserQFileSystemModelDialog()
    dialog.setDirRootPath(path)
    realPath = dialog.getRealPath(path)
    self.assertFalse("$" in realPath)
    newRealPath = os.path.join(realPath, "toto")
    newRealPathWithEnvVar = dialog.toPathWithEnvVar(newRealPath)
    self.assertTrue("$PACKAGESPY_ROOT_DIR" in newRealPathWithEnvVar)
    self.assertTrue("toto" in newRealPathWithEnvVar)

  def xtest_300(self):
    app = OnceQApplication()
    dialog = FISY.QFileSystemModelDialog()
    rootPath = "$PACKAGESPY_ROOT_DIR" 
    #os.path.realpath(os.path.join(testDir, '..', '..'))
    dialog.setDirRootPath(rootPath, filters=["*.py","*.c*","*.h*","*.C*"])
    dialog.show()
    self.launchTimer(app, dialog)

  def test_400(self):
    app = OnceQApplication()
    dialog = FISY.FileSystemModelViewerWidget()
    #path = "/volatile/wambeke/SAT4/prerequisites/uranie-3.6.0/FROM_nothing/share/uranie/html/index.html"
    #path = "helpWidgetUserQFileSystemModelDialog.html" TODO?
    #dialog.webView.load(QtCore.QUrl(path))
    #rootPath = os.path.join("$PACKAGESPY_ROOT_DIR", 'pythonAppliMatix')
    #rootPath = os.path.join("$URANIESYS", "share", "uranie")
    rootPath = os.path.join(testDir, "..", "..")
    #dialog.setDirRootPath(rootPath, filters=["*.py","*.c*","*.h*","*.C*"])
    dialog.setDirRootPath(rootPath, filters=["*.py"])
    dialog.show()
    self.launchTimer(app, dialog)

if __name__ == '__main__':
  verbose = False
  setTimer = False #user call and wait
  unittest.main()
  pass
