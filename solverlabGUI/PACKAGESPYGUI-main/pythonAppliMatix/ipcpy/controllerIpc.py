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
import pprint as PP
from PyQt5 import QtGui, QtCore, QtWidgets
import threading
from inspectpy.easyInspect import Inspect

import ipcpy.pipeFifo as PFIFO
from interpreterpy.interpreterPython import InterpreterPython as IPY
from salomepy.onceQApplication import OnceQApplication
from xyzpy.guiXyz.dialogXmlXyz import DialogXyz, QTextEditXyz

_aDir = os.path.realpath(os.path.join("/tmp", os.getenv("USER"), "testpipefifo"))
_aPipeName = os.path.realpath(os.path.join(_aDir, "pipefifo_test_001"))

verbose = False

########################################################
class ControllerIpc(QtCore.QObject):

  r'''
  controller of one thread for inter process communication (ipc)
  only read for command-messages and emit controllerOfModel pyqt signals
  - from external process 
  - to internal controllerOfModel instance

  # ONE:
  external process
    #as simple python code compliant, possibly multiline
    write(namedPipe, "IpcRefreshView()")

  # TWO:
  class ControllerIpc()
    def IpcRefreshViews(self):
      """refresh all views etc """
      self._controllerOfModel.refreshViewsSignal.emit()

  # THREE:
  class ControllerOfModel()
    refreshViewsSignal = QtCore.pyqtSignal()
    def __init__(self):
      self.refreshViewsSignal.connect(self.refreshModelViewsFromSignal)

    def refreshModelViewsFromSignal(self):
  '''

  '''
  note: 
    may be smarter for future
    http://sametmax.com/crossbar-le-futur-des-applications-web-python/
    https://flask-restful.readthedocs.io/en/latest/quickstart.html
  '''

  receiveSignal = QtCore.pyqtSignal(str)

  def __init__(self, *args, **kwargs):
    super(ControllerIpc, self).__init__(*args)
    self._className = self.__class__.__name__ #shortcut
    self._threadPipeIn = None
    self._controllerOfModel = None
    self.receiveSignal.connect(self.receiveMessage)
    self._results = []

  def __repr__(self):
    pipe = self._threadPipeIn
    if pipe == None:
      res = "%s(pipe=None)" % self.__class__.__name__
    else:
      res = "%s(pipe=%s, nbReceive=%i)" % \
            (self._className, pipe.getBaseNameFifo(), pipe.nbReceive)
    return res

  '''
  def __del__(self): #risky??? or useless???
    """http://enki-editor.org/2014/08/23/Pyqt_mem_mgmt.html"""
    print "%s__del__" % self.__class__.__name__
    self.stopThreadPipeIn() #as useless now
  '''

  def setControllerOfModel(self, controller):
    if self._controllerOfModel != None: #only once
      raise Exception("%s controllerOfModel set yet" % self._className)
    self._controllerOfModel = controller

  def setThreadPipeIn(self, nameFilePipeIn, Verbose=False):
    """create thread to receive, and start it"""
    if self._threadPipeIn != None:
      raise Exception("threadPipeIn set yet: %s" % self._threadPipeIn.getBaseNameFifo())
    self._threadPipeIn = PFIFO.ThreadPipeFifoReceiveBlocking(self, nameFilePipeIn)
    #print "TODO debug init ThreadPipeFifoBase"
    #self._threadPipeIn = PFIFO.ThreadPipeFifoBase()
    self._threadPipeIn.verbose = Verbose
    self._threadPipeIn.start()
    if Verbose: 
      print("setThreadPipeIn.start():\n  %s" % self._threadPipeIn.getNameFifo() )

  def stopThreadPipeIn(self, Verbose=False):
    if self._threadPipeIn != None:
      self._threadPipeIn.stop(Verbose=Verbose)

  def getThreadPipeIn(self):
    return self._threadPipeIn

  def receiveMessage(self, aMessage):
    mess = '%s receiveMessage "%s"' % (self._className, aMessage)
    #if self._threadPipeIn.verbose: print(mess)
    result = self.executeMessage(aMessage)
    self._results.append(mess + " result=%s" % result)

  def executeMessage(self, aMessage):
    ip = IPY( {"CIPC": self} ) #locals access to self.IpcSomething() methods
    theCode = str(aMessage)

    #for example...
    theExampleCode = r"""
# could be multiline code
CIPR__Repr__ = CIPC.__repr__()
bb = 'hello world %i' % (1+2)
"""

    ip.runcode(theCode) 
    if ip.isOk():
      if False: #self._threadPipeIn.verbose: 
        print("\n########### INFO: ControllerIpc:\n%s\n" % ip.getStrResume(theCode))
      return True
    else:
      return self.createEditorSource(theCode, ip.getStdErr(), parent=None)

    """else: #single message
      if self._threadPipeIn.verbose:
        print("\n########### ERROR: ControllerIpc:\n%s\n" % ip.getStrResume(theCode))
      return False"""

  def createEditorSource(self, theCode, error, parent=None):
    if parent == None: pparent=self._controllerOfModel.getDesktop()
    aDialog = DialogXyz(parent=pparent) # open in center desktop
    aWidget = QTextEditXyz()
    aWidgetErr = QTextEditXyz()
    aDialog.setUpWidgetLayout([aWidget, aWidgetErr])
    #aWidgetErr.hide()
    aWidget.setValue(theCode)
    aWidgetErr.setValue(error)
    aDialog.setMinimumSize(500, 450)
    aDialog.setWindowTitle("uranieGui IPC problem")
    while True:
      aDialog.exec_()
      if aDialog.choice == "Cancel": return True

      newSource = aWidget.getValue()
      if verbose: 
        print("createEditorSource: new source:'\n%s'" % newSource)
      
      #test newSource is correct
      ip = IPY( {"CIPC": self} ) #locals access to self.IpcSomething() methodsip = 
      ip.runcode(newSource)
      if ip.isOk(): break
      #print "Ooops:\n",ip.getStdErr()
      aWidgetErr.setValue(ip.getStdErr())
      aWidgetErr.show()
      # TODO fix why??? if not QMessageBox, next aDialog.exec_() is not in center desktop
      QtWidgets.QMessageBox.warning(self.getDesktop(), "error in source", ip.getStdErr())

    return True


  def printResults(self):
    """only for multithread unittest"""
    print("\n%s results:\n%s\n"  % (self.__class__.__name__, PP.pformat(self._results)))

  def printApiIpcDoc(self):
    """idem IpcprintApi"""
    self.IpcPrintApi()

  def getApiIpcDoc(self):
    mess = Inspect(self).get_prefixedMethodsDoc("Ipc")
    return mess

  def getApiDoc(self):
    mess = Inspect(self).get_prefixedMethodsDoc("")
    return mess

  ######################################################
  # elementaries IpcMethods
  # user have inherited class for his particular stuff
  # see class ControllerIpcUra, for example
  #####################################################
  def IpcRefreshViews(self):
    """
    refresh GUI TreeView
    """
    self._controllerOfModel.refreshViewsSignal.emit()
    
  def IpcPrintApi(self):
    """
    print authorized ControllerIpc API (and so inherited classes)
    known as all methods (with 'Ipc' prefix) 
    for external process sender ipc use
    """
    mess = "\n#### %s: Ipc API:\n%s" % (self.__class__.__name__, self.getApiIpcDoc())
    print(mess)
    

########################################################
class ControllerOfModel(QtCore.QObject):
  """
  Class as ControllerXyz (for test and user example)
  simple for demo unittest of ControllerIpc
  """
  refreshViewsSignal = QtCore.pyqtSignal()
  doSomethingSignal = QtCore.pyqtSignal(dict) 
  # ...
  
  def __init__(self):
    super(ControllerOfModel, self).__init__()
    self.refreshViewsSignal.connect(self.refreshModelViewsFromSignal)
    self.doSomethingSignal.connect(self.doSomethingFromSignal)
    self._results = []
    self.verbose = False

  def refreshModelViewsFromSignal(self):
    mess = "%s execute refreshModelViewsFromSignal" % self.__class__.__name__
    if self.verbose: print(mess)
    self._results.append(mess)

  def doSomethingFromSignal(self, param1="param1default", param2="param2default"):
    """
    doSomethingFromSignal obviously do something
    - param1 = ...
    - param2 = ...
    """
    mess = "%s execute doSomethingFromSignal param1=% param2=%s" % (self.__class__.__name__, param1, param2)
    if self.verbose: print(mess)
    self._results.append(mess)

  def printResults(self):
    """only for multithread unittest"""
    print("\n%s results:\n%s\n"  % (self.__class__.__name__, PP.pformat(self._results)))



