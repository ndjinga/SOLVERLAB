#!/usr/bin/env python
# -*- coding: utf-8 -*-

# %% LICENSE_SALOME_CEA_BEGIN
# see PACKAGESPY/LICENCE file
# %% LICENSE_END


import unittest
import sys
import subprocess as SP
from salomepy.onceQApplication import OnceQApplication
from salomepy.qTextEditForLog import QTextEditForLog
from salomepy.threadWorkerForWidgetEdit import ThreadWorkerForWidgetEdit

from PyQt5 import QtCore, QtGui, QtWidgets

verbose = False
setTimer = True
timers = []
deltaTime = 6000
withShow = True

def aTestQTextEditForLog():
  """only for detailled function through ThreadWorker, see elementary test_010"""
  setTimer=True #False #user call and wait
  app = OnceQApplication()
  cmd = ["python", "-u", "-c", 
          "print 'An stdout hello word in black from python Popen!' ;" +
          "raise Exception('An stderr Error example in red From python!')" ]
  proc1 = SP.Popen(cmd, stdout=SP.PIPE, stderr=SP.PIPE)
  edit = QTextEditForLog()
  edit.show()
  edit.stdout_worker1 = ThreadWorkerForWidgetEdit(proc1.stdout, edit, "Black")
  edit.stderr_worker1 = ThreadWorkerForWidgetEdit(proc1.stderr, edit, "Red")
  edit.stdout_worker1.start()
  edit.stderr_worker1.start()
  cmd = "echo 'some text in blue...'"
  cmd += "; sleep 2"
  cmd += "; echo 'll ll  llllll  l monospaced...'"
  cmd += "; echo 'mm mm2bmmmmmm  m or not?...'"
  cmd += "; ls -alt *"
  edit.proc2 = SP.Popen(cmd, shell=True, stdout=SP.PIPE, stderr=SP.PIPE)
  edit.stdout_worker2 = ThreadWorkerForWidgetEdit(edit.proc2.stdout, edit, "Blue")
  edit.stderr_worker2 = ThreadWorkerForWidgetEdit(edit.proc2.stderr, edit, "Magenta")
  edit.stdout_worker2.start()
  edit.stderr_worker2.start()
  for i in edit.getCurrentFontFamilies():
    #first arrived, first shown
    edit.insertLine(i)
  if verbose: 
    print("QTextEditForLog Main is finished, but not thread popened through QTextEditForLog instance!")
  return edit

class TestCase(unittest.TestCase):

  def launchTimer(self, wid):
    if setTimer: 
      self.app = OnceQApplication()
      timer = QtCore.QTimer();
      timer.timeout.connect(wid.close)
      timer.start(deltaTime)
      timers.append(timer)

  def xtest_000(self):
    app = OnceQApplication()
    edit = aTestQTextEditForLog()
    self.launchTimer(edit)
    app.exec_()
  
  def test_010(self):
    app = OnceQApplication()
    cmd = "echo 'some text... and sleep 1'"
    cmd += "; sleep 1" #have not to stop threads on closeEvent
    cmd += "; ls -alt *"
    cmd += "; someExampleError..."
    edit = QTextEditForLog()
    edit.launchIntoPopen(cmd)
    edit.show()
    self.launchTimer(edit)
    app.exec_()
 
  def test_020(self):
    app = OnceQApplication()
    cmd = "echo 'some text... and sleep 12'"
    cmd += "; sleep 12" #have to stop threads on closeEvent
    cmd += "; ls -alt *"
    cmd += "; someExampleError..."
    edit1 = QTextEditForLog()
    edit1.launchIntoPopen(cmd)
    edit1.show()
    edit2 = QTextEditForLog()
    edit2.launchIntoPopen(cmd)
    edit2.show()
    self.launchTimer(edit1)
    self.launchTimer(edit2)
    app.exec_()

  def test_030(self):
    app = OnceQApplication()
    cmd = "echo 'some text... and sleep %i'"
    cmd += "; sleep %i" #have to stop threads on closeEvent
    cmd += "; someExampleError_%i..."
    edit1 = QTextEditForLog()
    for i in range(0, 10):
      edit1.launchIntoPopen(cmd % (i, i, i))
    edit1.show()
    self.launchTimer(edit1)
    app.exec_()
    
  def test_999(self):
    timers = [] #garbage collecting
    
if __name__ == '__main__':
  setTimer=True #True #False #user call and wait
  unittest.main()
  

