#!/usr/bin/env python
# -*- coding: utf-8 -*-

# %% LICENSE_SALOME_CEA_BEGIN
# see PACKAGESPY/LICENCE file
# %% LICENSE_END

import os
import unittest

from PyQt5 import QtCore, QtGui, QtWidgets
from salomepy.onceQApplication import OnceQApplication
from matplotlibpy.matplotlibWindowToolbar import MatplotlibWindowToolbar
from pylab import *
import salomepy.utilsWorkdir as UTW

testDir = UTW.getWorkdirDefault("TESTS")  

verbose = False
setTimer = True
timers = []
deltaTime = 500
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

  def test_010(self):
    app = OnceQApplication()
    fen = MatplotlibWindowToolbar()
    fig = fen.getFigure()
    fig.clear()
    ax = fig.add_subplot(1, 1, 1) 
    ax.plot(self.x,self.y1)
    ax.plot(self.x,self.y2)
    ax.set_title(self.title)
    ax.set_xlabel(self.xlabel)
    ax.set_ylabel(self.ylabel)
    fen.show()
    self.launchTimer(fen)
    app.exec_()

  def test_020(self):
    app = OnceQApplication()
    fen = MatplotlibWindowToolbar()
    fen.set_curve(self.x, self.y1, self.title, self.subtitle, self.xlabel, self.ylabel)
    fen.show()
    self.launchTimer(fen)
    app.exec_()

  def test_030(self):
    app = OnceQApplication()
    fen = MatplotlibWindowToolbar()
    fen.set_2_curves(self.xy12, self.title, self.subtitle, self.xy12labels)
    fen.show()
    self.launchTimer(fen)
    app.exec_()

  def test_040(self):
    app = OnceQApplication()
    fen = MatplotlibWindowToolbar()
    fen.set_n_curves(self.xy13, self.title, self.subtitle, self.xy13labels, title_y="y_label_unity")
    fen.show()
    self.launchTimer(fen)
    app.exec_()
    
  def test_045(self):
    app = OnceQApplication()
    fen = MatplotlibWindowToolbar()
    afile = os.path.join(testDir, "maplotlibTest.csv")
    s = ""
    for i in range(10):
      ii = float(i)
      s += "%i %.5e %.4f %s\n" % (ii*2+1, (ii+10)/2, ii**1.5, str((11-ii)**.5))
    with open(afile, "w") as F: F.write(s)
    fen.set_n_plots_from_files([afile])
    fen.show()
    self.launchTimer(fen)
    app.exec_()
    
  def test_050(self):
    app = OnceQApplication()
    fen = MatplotlibWindowToolbar()
    fen.show()
    fen.allTests()
    self.launchTimer(fen)
    app.exec_()

    
  def test_999(self):
    timers = [] #garbage collecting
    
if __name__ == '__main__':
  setTimer=False #True #False #user call and wait
  deltaTime = 2000
  unittest.main()
  

