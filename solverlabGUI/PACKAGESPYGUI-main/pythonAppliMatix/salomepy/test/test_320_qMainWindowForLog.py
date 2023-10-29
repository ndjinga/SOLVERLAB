#!/usr/bin/env python
# -*- coding: utf-8 -*-

# %% LICENSE_SALOME_CEA_BEGIN
# see PACKAGESPY/LICENCE file
# %% LICENSE_END


import unittest

from PyQt5 import QtCore, QtGui, QtWidgets

from salomepy.onceQApplication import OnceQApplication
from salomepy.qMainWindowForLog import QMainWindowForLog


verbose = False
setTimer = True
timers = []
deltaTime = 6000
withShow = True

class TestCase(unittest.TestCase):

  def launchTimer(self, wid):
    if setTimer: 
      self.app = OnceQApplication()
      timer = QtCore.QTimer();
      timer.timeout.connect(wid.close)
      timer.start(deltaTime)
      timers.append(timer)

  def test_010(self):
    app = OnceQApplication()
    cmd = "echo 'some text... and sleep 1'"
    cmd += "; sleep 1" #have not to stop threads on closeEvent
    cmd += "; ls -alt *"
    cmd += "; someExampleError..."
    fen = QMainWindowForLog()
    fen.launchCmdIntoPopen(cmd)
    fen.show()
    self.launchTimer(fen)
    app.exec_()

  def test_020(self):
    app = OnceQApplication()
    cmd = "echo 'some text... and sleep 1'"
    cmd += "; sleep 1" #have not to stop threads on closeEvent
    cmd += "; ls -alt *"
    cmd += "; someExampleError..."
    fen = QMainWindowForLog( withDocks=False )
    fen.launchCmdIntoPopen(cmd)
    fen.show()
    self.launchTimer(fen)
    app.exec_()
    
  def test_999(self):
    timers = [] #garbage collecting
    
if __name__ == '__main__':
  setTimer=True #True #False #user call and wait
  unittest.main()
  

