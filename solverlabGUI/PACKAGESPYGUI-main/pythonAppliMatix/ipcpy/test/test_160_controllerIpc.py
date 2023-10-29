#!/usr/bin/env python
# -*- coding: utf-8 -*-

# %% LICENSE_SALOME_CEA_BEGIN
# see PACKAGESPY/LICENCE file
# %% LICENSE_END

r'''
All methods named ControllerIpc.IpcSomething are allowed 
as external process commands 
to be executed by internal associated controllerOfModel methods
called throught pyqt controllerOfModel signals
'''

import os
import stat
import time
import pprint as PP
from PyQt5 import QtCore, QtWidgets
import unittest
import threading

from ipcpy.controllerIpc import ControllerOfModel, ControllerIpc
from inspectpy.easyInspect import Inspect

import ipcpy.pipeFifo as PFIFO
from salomepy.onceQApplication import OnceQApplication

_aDir = os.path.realpath(os.path.join("/tmp", os.getenv("USER"), "testpipefifo"))
_aPipeName = os.path.realpath(os.path.join(_aDir, "pipefifo_test_001"))

verbose = False
deltaTime = 10000 # milliseconds max delay for threaded tests


########################################################
class TestCase(unittest.TestCase):  
  
  def toHex(self, aStr):
    aHex = ":".join("{:02x}".format(ord(c)) for c in aStr)
    return aHex
  
  def printHex(self, aStr):
    print('"%s" -> \n%s' % (aStr, self.toHex(aStr)))

  def launchApp(self, theTest, widget=None):
    """
    launch blocking app.exec_()
    MUST be called in the main thread
    for pyqt signals events
    launch threaded theTest before app.exec_()
    """ 
    aThread = threading.Thread(target=theTest)
    aThread.daemon = True
    app = OnceQApplication()
    if widget == None:
      wid = QtWidgets.QLabel("from %s\n\nClose me (x) for exit Qt event loop !\n" % os.path.realpath(__file__))
      # or QtWidgets.QMainWindow()
    else:
      wid = widget
    wid.show()
    timer = QtCore.QTimer()
    timer.timeout.connect(wid.close)
    timer.start(deltaTime)
    aThread.start()
    if verbose: print("!!! begin QApplication loop !!!")
    app.exec_()  # code blocks over here if wid no closed
    if verbose: print("!!! end QApplication loop !!!")
    return aThread
    
  def test_000(self):
    try:
      os.remove(_aPipeName)
    except:
      pass
    self.assertFalse(os.path.exists(_aPipeName))

  def test_005(self):
    # test of Inspect
    if verbose: print("\n**** test_005 ****\n")
    a = ControllerOfModel()
    b = ControllerIpc()
    stra = Inspect(a).get_pydoc_render_doc()[:200]
    #self.printHex(stra)
    self.assertTrue( r"ControllerOfModel(PyQt5.QtCore.QObject)" in stra )

    strb = Inspect(b).get_pydoc_render_doc()[:200]
    #self.printHex(strb)
    self.assertTrue( r"class ControllerIpc(PyQt5.QtCore.QObject)" in strb )

    self.assertTrue( \
       'doSomethingFromSignal(self, param1="param1default", param2="param2default")' in \
       Inspect(a).get_inspect_signature(a.doSomethingFromSignal) )
    #self.assertTrue( \
    #   'IpcLoadDataFile(self, aFile=None)' in \
    #    Inspect(b).get_inspect_signature(b.IpcLoadDataFile) )

    apidoc = b.getApiIpcDoc()
    if verbose: print("\n#### ControllerIpc.getApiIpcDoc():\n%s" % apidoc)
    #self.assertTrue( 'METHOD: IpcLoadDataFile(self, aFile=None)' in apidoc)
    self.assertTrue( 'METHOD: IpcRefreshViews(self)' in apidoc)

  def test_010(self):
    if verbose: print("\n**** test_010 ****\n")
    a = PFIFO.PipeFifo()
    a.setNameFifo(_aPipeName)
    self.assertTrue(a.isPipeFifo())

  def test_020(self):
    if verbose: print("\n**** test_020 ****\n")
    cipc = ControllerIpc()
    # TODO more...

  def test_025(self):
    if verbose: print("\n**** test_025 ****\n")
    com = ControllerOfModel()
    # TODO more...

  def threaded_test_030(self):
    if verbose: print("\n**** begin threaded_test_030 ****\n")
    thSend = PFIFO.ThreadPipeFifoSend(None, _aPipeName)
    thSend.start()
     # main and PipeIn and launchApp.aThread and thSend
    self.assertEqual(len(PFIFO.getListOfThreadsAlive()), 4)
    if verbose: thSend.sendMessage("CIPC.IpcPrintApi()")
    time.sleep(2)
    thSend.sendMessage("CIPC.IpcRefreshViews()")
    #time.sleep(2)
    #thSend.sendMessage("CIPC.IpcLoadDataFile(aFile='/tmp/aDataFile')")
    thSend.stop()
    time.sleep(2)
    if verbose: print("\n**** end threaded_test_030 ****\n")
    
  def threaded_test_040(self):
    print("\n**** threaded_test_040 ****\n")
    time.sleep(1)
    self.assertEqual(len(PFIFO.getListOfThreadsAlive()), 2)
    
  def test_900(self):
    # threaded_tests
    if verbose: print("\n**** test_900 ****\n")
    cipc = ControllerIpc()
    cipc.verbose = verbose
    cipc.setParent(OnceQApplication())
    com = ControllerOfModel()
    cipc.setControllerOfModel(com)
    cipc.setThreadPipeIn(_aPipeName)
    # main and PipeIn
    self.assertEqual(len(PFIFO.getListOfThreadsAlive()), 2)
    
    self.launchApp(self.threaded_test_030)
    if verbose: 
      cipc.printResults()  
      com.printResults()
    cipc.stopThreadPipeIn()
    time.sleep(2)

    #self.launchApp(self.threaded_test_040)    

  def xtest_040(self):
    if verbose: print("\n**** test_040 ****\n")
    app = OnceQApplication() #created in the main() thread
    thApp = self.launchApp()
    time.sleep(1)
    self.assertEqual(len(PFIFO.getListOfThreadsAlive()), 1)


if __name__ == '__main__':
  verbose = False #True #False
  PFIFO.verbose = verbose # True #
  setTimer = not verbose #user call and wait
  unittest.main()
  pass
