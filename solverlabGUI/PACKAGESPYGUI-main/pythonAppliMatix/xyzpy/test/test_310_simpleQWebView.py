#!/usr/bin/env python
# -*- coding: utf-8 -*-

# %% LICENSE_SALOME_CEA_BEGIN
# see PACKAGESPY/LICENCE file
# %% LICENSE_END



import os
import sys
import unittest

from PyQt5 import QtCore, QtGui, QtWidgets
from salomepy.onceQApplication import OnceQApplication
from xyzpy.guiXyz.FileSystemXyz import BrowserWebView as Browser

setTimer = True
deltaTime = 1000
verbose = False
verboseView = False
verboseEvent = False

testDir = os.path.split(os.path.realpath(__file__))[0]

# print testDir

############################################
class TestCase(unittest.TestCase):
  
  def launchTimer(self, app, wid):
    if setTimer: 
      timer = QtCore.QTimer();
      timer.timeout.connect(wid.close)
      timer.start(deltaTime)
    app.exec_()

  def test_400(self):
    import salomepy.utilsWorkdir as UTW
    app = OnceQApplication()
    wid = Browser()
    wid.resize(900, 500)
    path = UTW.getHtmlDoc()
    # generate (No document for ...' if not os.path.isfile(path)
    if path is not None: wid.load(QtCore.QUrl(path))
    wid.show()
    self.launchTimer(app, wid)

if __name__ == '__main__':
  verbose = False
  setTimer = False #user call and wait
  unittest.main()
  pass
