#!/usr/bin/env python
# -*- coding: utf-8 -*-

# %% LICENSE_SALOME_CEA_BEGIN
# see PACKAGESPY/LICENCE file
# %% LICENSE_END


import unittest
import sys
from PyQt5 import QtCore, QtGui, QtWidgets
from salomepy.onceQApplication import OnceQApplication
from salomepy.qTabMultipleTextCentral import QTabMultipleTextCentral
from time import sleep

verbose = False
setTimer = True
deltaTime = 3000
withShow = True

class TestCase(unittest.TestCase):

  def launchTimer(self, wid):
    self.app = OnceQApplication()
    self.timer = QtCore.QTimer();
    self.timer.timeout.connect(wid.close)
    if setTimer: self.timer.start(deltaTime)
    self.app.exec_()

  def test_010(self):
    app = OnceQApplication()
    fen = QTabMultipleTextCentral()
    fen.logSalomeWidget.launchIntoPopen("pwd ; echo sleep 2 ; sleep 2 ; ls -alt *py ; sommmeErrorExampleInRed")
    fen.showLogSalomeWidget()
    fen.show()
    self.launchTimer(fen)
    
  def test_020(self):
    app = OnceQApplication()
    fen = QTabMultipleTextCentral()
    #closeEvent and stop threads
    fen.logSalomeWidget.launchIntoPopen("pwd ; echo sleep 10 and closeEvent... ; sleep 10 ; ls -alt *py ; sommmeErrorExampleInRed")
    fen.showLogSalomeWidget()
    fen.show()
    self.launchTimer(fen)
    
  
if __name__ == '__main__':
  setTimer=True #False #user call and wait
  unittest.main()
  pass

  

