#!/usr/bin/env python
# -*- coding: utf-8 -*-

# %% LICENSE_SALOME_CEA_BEGIN
# see PACKAGESPY/LICENCE file
# %% LICENSE_END


import unittest

from PyQt5 import QtCore, QtGui, QtWidgets
from salomepy.onceQApplication import OnceQApplication
from salomepy.browserPythonClass import *
from matplotlibpy.matplotlibWindowToolbar import MatplotlibWindowToolbar
from pylab import *


verbose = False
setTimer = True
timers = []
deltaTime = 1000
withShow = True

class TestCase(unittest.TestCase):
  x = arange(0.0, 2.0, 0.01)
  y1= sin(2*pi*x)
  y2= .5*sin(3*pi*x)
  y3= .2*sin(6*pi*x)
  
  title = "a_title"
  subtitle = "a_sub_title"
  xlabel = "x_label"
  ylabel = "y_label"
  y1label = "y1_label"
  y2label = "y2_label"
  y3label = "y3_label"
  
  xy12 = [x, y1, y2]
  xy12labels = [xlabel, y1label, y2label ]
  xy13 = [x, y1, y2, y3]
  xy13labels = [xlabel, y1label, y2label , y3label ]

  def launchTimer(self, wid):
    if setTimer: 
      self.app = OnceQApplication()
      timer = QtCore.QTimer();
      timer.timeout.connect(wid.close)
      timer.start(deltaTime)
      timers.append(timer)

  def xtest_010(self):
    a = 1
    printDir(a)

  def test_020(self):
    app = OnceQApplication()
    a = 1
    #for i in dir(a): print "dir:", i
    res = BrowserPythonClass(a)
    if verbose: print(res.toStrXml())
    fen = res.display()
    self.launchTimer(fen)
    app.exec_()

  def test_022(self):
    app = OnceQApplication()
    a = [1, 2] 
    res = BrowserPythonClass(a)
    if verbose: print(res.toStrXml())
    fen = res.display()
    self.launchTimer(fen)
    app.exec_()

  def test_024(self):
    app = OnceQApplication()
    a = (1, 2)
    res = BrowserPythonClass(a)
    if verbose: print(res.toStrXml())
    fen = res.display()
    self.launchTimer(fen)
    app.exec_()

  def test_026(self):
    app = OnceQApplication()
    a = {"key1": "aKey1Value", 22: "aKey22Value", 33.: "aKey22Value"}
    res = BrowserPythonClass(a)
    if verbose: print(res.toStrXml())
    fen = res.display()
    self.launchTimer(fen)
    app.exec_()

  def test_030(self):
    app = OnceQApplication()
    a = 1
    res = BrowserPythonClass(a)
    res = BrowserPythonClass(res, levelMax=4)
    fen = res.display()
    self.launchTimer(fen)
    app.exec_()

  def test_040(self):
    app = OnceQApplication()
    a = 1
    res = BrowserPythonClass(a)
    res.anExempleDict = {"key1": "aKey1Value", 22: "aKey22Value"}
    res.anExempleList = ["hello", 123, [1, 2, 3.]]
    res = BrowserPythonClass(res, levelMax=2)
    fen = res.display()
    self.launchTimer(fen)
    app.exec_()

  def test_050(self):
    #avoid levelMax=4 too long?
    app = OnceQApplication()
    fen = MatplotlibWindowToolbar()
    fen.set_curve(self.x, self.y1, self.title, self.subtitle, self.xlabel, self.ylabel)
    res = BrowserPythonClass(fen.getFigure()._axstack._elements, levelMax=1)
    fen.show()
    fen2 = res.display()
    self.launchTimer(fen)
    self.launchTimer(fen2)
    app.exec_()
    
    
if __name__ == '__main__':
  setTimer=True #False #user call and wait
  unittest.main()
  

