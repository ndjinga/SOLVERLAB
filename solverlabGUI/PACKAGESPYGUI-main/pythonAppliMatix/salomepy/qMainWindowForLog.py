#!/usr/bin/env python
# -*- coding: utf-8 -*-

# %% LICENSE_SALOME_CEA_BEGIN
# see PACKAGESPY/LICENCE file
# %% LICENSE_END


import os
import sys
import shutil
import fnmatch
import glob
from time import sleep
from datetime import datetime

from PyQt5 import QtGui, QtCore, QtWidgets
from salomepy.strEvent import *
from salomepy.qTabMultipleTextCentral import QTabMultipleTextCentral
import salomepy.iconsUser as IUSR
from xyzpy.guiXyz.treeXmlXyz import TreeXmlXyz
import xyzpy.loggingXyz as LOG

logger = LOG.getLogger()

verbose = False
verboseEvent = False


class QMainWindowForLog(QtWidgets.QMainWindow):  
  # Define a new signal called 'launchOnArgs'
  launchOnArgs = QtCore.pyqtSignal()
  index = [0]

  def __init__(self, *args, **kwargs):
    super(QMainWindowForLog, self).__init__(*args)
    self.args = args
    self.kwargs = kwargs
    if verbose: print("QMainWindowForLog args %s, kwargs %s" % (str(args), str(kwargs)))
    self.setObjectName("QMainWindowForLog"+str(self.index))
    #self.setAttribute(QtCore.Qt.WA_DeleteOnClose, False)
    self.index[0] += 1 #unambigous objectName
    self.setWindowTitle(self.objectName())
    self.setWindowModality(QtCore.Qt.NonModal)
    self.prefixShortcut = "Ctrl+"
    self.treeView = None
    self.__addCentral()
    self.__addDocks()
    self.__createActions()
    self.__addToolBars()
    self.statusBar().showMessage('Ready')
    self.resize(900, 500)
    self.salomeStudy = None
    #self.myParent = myParent #may be not a QWidget... to reach methods test of myParent
    self.currentDir = None
    self.currentDirResultats = None
    self.currentMed = None
    self.currentSauv = None
    self.currentDgibi = None
    self.doneCases=[]
    self.launchOnArgs.connect(self.theLaunchOnArgs)
    self.treeView.quickEditFiles.connect(self.centralWidget().quickEditFiles)
    #self.centralWidget().logCastemWidget.finishPopen.connect(self.treeView.refreshModel)
    self.initializeWorkdirDone = False
    self.initializeWorkdir()
    #self.controller.connectTreeViewAndTableEdit(self.treeView, self.centralWidget().csvEditWidget)
    self._controller = None
    if "withDocks" in self.kwargs: 
      withDocks = bool(self.kwargs["withDocks"])
    else:
      withDocks = True
    if withDocks == True:
      for dock in self.docks: dock.show()
    else:
      for dock in self.docks: dock.hide()
  
  """
  #warning: is not usable if salome, only catch hide event, not close
  def close(self):
    print("QMainWindowForLog %s close" % self.objectName())
    return super(QMainWindowForLog, self).close()

  def closeEvent(self, event):
    print("QMainWindowForLog %s closeEvent" % self.objectName())
    #event.ignore()
    return super(QMainWindowForLog, self).closeEvent(event)

  def event(self, event):
    print("QMainWindowForLog %s event" % (self.objectName(),strEvent(event)))
    return super(QMainWindowForLog, self).event(event)

  def __del__(self):
    print("QMainWindowForLog %s __del__" % self.objectName())
    return super(QMainWindowForLog, self).__del__()
  """

  def setElementaryLog(self):
    """hide tree and no log widget"""
    for dock in self.docks: dock.hide()
    ctabs = self.centralWidget()
    ctabs.hideOtherTexteditWidget()
    ctabs.hideXmlWidget()
  
  def setController(self, controller):
    """really could be None if no use in view without MVC pattern"""
    if self._controller == None:
      self._controller = controller
      return
    raise Exception("QMainWindowForLog.setController done yet for %s as %s" % (self.objectName(), self.getController().objectName()))
    
  def getController(self):
    """to get (for example) sendRequest method of controller api"""
    return self._controller

  def receiveRequestToView(self, strXmlRequest):
    if verboseEvent: 
      print("%s %s receiveRequest virtual" % (self.__class.__name__, self.objectName()))
    return True

  def getIncrement(self):
    """
    get an once shot useful time sorted str, based on datetime.now()
    like numeric YearMonthDayHour...
    """
    now = datetime.now()
    #return now.strftime("%Y-%m-%d %H:%M")
    #return now.strftime("%m%d%H%M%S")
    return now.strftime("%H%M%S")

  def __addCentral(self):
    try:
      central = self.kwargs["centralWidget"]
    except:
      central = QTabMultipleTextCentral()

    self.setCentralWidget(central)
    #self.centralWidget().resize(self.centralWidget().size())
    self.centralWidget().show()

  def getTabByName(self, aName):
    """
    name as attribute name or name as tab label
    wid.tabName = "Explore Dir"      #tab label
    wid.attName = "exploreDirWidget" #attribute name
    """
    return self.centralWidget().getTabByName(aName)
    
  def getTabs(self):
    return self.centralWidget().getTabs()

  def showTab(self, index):
    self.centralWidget().showTab(index)

  def showTabByName(self, name):
    self.centralWidget().showTabByName(name)

  def __addDocks(self):
    self.docks = []
    dock = QtWidgets.QDockWidget("TreeForLog", self)
    self.treeView=TreeXmlXyz()
    logger.debug("create treeView")
    dock.setWidget(self.treeView)
    dock.setAllowedAreas(QtCore.Qt.LeftDockWidgetArea)
    self.docks.append(dock)
    for dock in self.docks:
      self.addDockWidget(QtCore.Qt.LeftDockWidgetArea, dock)
    
  def __addToolBars(self):
    self.toolBars = []
    tb = self.addToolBar("EditLog")
    for action in self.actions:
      tb.addAction(action)
    self.toolBars.append(tb)
  
  def __createActions(self):
    self.actions = []
    #self.actions.append(self.__createAction( "Dgibi", "D", "Show file dgibi", self.showDgibiAction))
    #self.actions.append(self.__createAction( "Castem", "C", "Launch Castem calculus", self.castemLaunchAction, "castem"))
    #self.actions.append(self.__createAction( "Paraview", "S", "Launch Paraview pattern", self.showResultsAction, "paraview"))
    self.actions.append(self.__createAction( "AllTestLauncher", "T", "Launch all python tests", self.testAllAction, "tests"))
    #self.actions.append(self.__createAction( "Cassis", "S", "Launch Cassis-Salome cession", self.getSalomeStudy))

  def __createAction(self, Name, Shortcut, ToolTip, Call, Icon=None):
    action = QtWidgets.QAction(Name, self)
    if Shortcut!=None: action.setShortcut(self.prefixShortcut+Shortcut)
    action.setToolTip(ToolTip)
    if Icon!=None:
      action.setIcon(IUSR.getIconFromName(Icon))
    action.triggered.connect(Call)
    return action

  def makedirs(self, namedir):
    if os.path.exists(namedir):
      dirbak=namedir+".bak"
      if os.path.exists(dirbak): shutil.rmtree(dirbak)
      os.rename(namedir,dirbak)
      os.listdir(dirbak) #sert seulement a mettre a jour le systeme de fichier sur certaines machines
    os.makedirs(namedir)
  
  def initializeWorkdir(self):
    """
    initialize user logs directory $WORKDIR4LOG if not existing
    """
    workDir=os.getenv("WORKDIR4LOG")
    if workDir == None: 
      workDir=os.path.join("/tmp",  os.getenv('USERNAME'))
    workDir = os.path.realpath(workDir)
    os.environ["WORKDIR4LOG"] = workDir
    if not os.path.exists(workDir):
      print('create inexisting $WORKDIR4LOG',workDir)
      logger.warning('create inexisting $WORKDIR4LOG ' + workDir)
      os.makedirs(workDir)
    logger.debug("QMainWindowForLog logs files in $WORKDIR4LOG '%s'" % workDir)
    #self.treeView.refresh.emit()
    
  def med2sauv(self, file_in, file_out=None):
    import medutilities as MU #needs import _MEDLoader salome libmedloader.so so not in main
    #if MU == None:
    #  logger.error("can't import medutilities, no med2sauv conversion")
    #  return None
    if file_out==None:
      dire, base=os.path.split(file_in) 
      file_out = os.path.join(dire, os.path.splitext(base)[0] + ".sauv")
    if not self.isSameFile(file_in, file_out):
      MU.convert(file_in, "MED", "GIBI", file_out=file_out)
      return file_out
    else:
      logger.error("med2sauv: not done: files are same: "+file_in+" -> "+file_out)
      return None
    
  def sauv2med(self, file_in, file_out=None):
    import medutilities as MU #needs import _MEDLoader salome libmedloader.so so not in main
    #if MU == None:
    #  logger.error("can't import medutilities, no sauv2med conversion")
    #  return None
    if file_out==None:
      dire, base=os.path.split(file_in) 
      file_out = os.path.join(dire, os.path.splitext(base)[0] + ".med")
    if not self.isSameFile(file_in, file_out):
      MU.convert(file_in, "GIBI", "MED", file_out=file_out)
      return file_out
    else:
      logger.error("sauv2med: not done: files are same: "+file_in+" -> "+file_out)
      return None

  def isSameFile(self, file1, file2):
    """follow relative paths and symbolic links"""
    rfile1=os.path.realpath(file1)
    rfile2=os.path.realpath(file2)
    return (rfile1==rfile2)

  def theLaunchOnArgs(self):
    """launch from command line"""
    #logger.debug("args.parsed_options: " + str(self.args.parser.values)) #,self.args.parse_args()
    aCase=PVU.PvCase(self.args.getAllOptionsAsDict())
    logger.info("theLaunchOnArgs " + aCase.action)
    if aCase.action=="cmdLaunch":
      self.cmdLaunchAction(aCase)

  def launchCmdIntoPopen(self, cmd, clearBefore=True, verbose=False, signal=None, cwd=None):
    if verbose: print("QMainWindowForLog.launchCmdIntoPopen command '%s'" % cmd)
    self.centralWidget().showLogCMDWidget()
    if clearBefore: self.centralWidget().logCMDWidget.clear()
    self.centralWidget().logCMDWidget.launchIntoPopen(cmd, signal=signal, cwd=cwd)

  def insertText(self, text, clearBefore=False):
    if clearBefore: self.centralWidget().logCMDWidget.clear()
    self.centralWidget().logCMDWidget.insertText(text)
  
  def testAllAction(self):
    logger.info("testAllAction")
    """ 
    import distenepy.utilsDisteneHelp as UDIS
    data = UDIS.hybridHelp()
    self.treeView.setFromXml(data.toXml())
    self.centralWidget().xmlWidget.insertText(data.toStrXml())
    """
    #unittests, for big example, AllTestLauncher.sh have to be in $PATH
    cmd = "cd ${PACKAGESPY_ROOT_DIR}; pwd ; AllTestLauncher.sh"
    self.launchCmdIntoPopen(cmd)
    
  def launchSalome(self):
    cmd = 'python -u -c "from salomeRunner import SalomeRunner ; runner = SalomeRunner(configFileNames=None) ; runner.go([])" &'
    self.launchSalomeIntoPopen(cmd)
    sleep(1) #with a background "&" the first sleep is killed ...mystery...
    self.salomeStudy = None

  def getSalomeStudy(self):
    """store result in self.salomeStudy"""
    if self.salomeStudy!=None:
      #if salome crash or GUI user exit
      try:
        logger.info("GetSalomeStudy: (done yet) got study: %s Id: %s",self.salomeStudy._get_Name(),self.salomeStudy._get_StudyId())
        return
      except:
        self.salomeStudy==None
    self.omniOrbConfigFile = os.getenv('OMNIORB_CONFIG')
    if self.omniOrbConfigFile==None:
      logger.error("I can't get a salome cession because OMNIORB_CONFIG environment variable is None.")
      self.salomeStudy = None
      return
    if not os.path.isfile(self.omniOrbConfigFile):
      logger.warning("I can't see a salome cession: no file " + self.omniOrbConfigFile)
      logger.warning("I try to launch a salome cession... take a time...")
      self.launchSalome()
    import salome
    for i in range(0,3): #try twice
      sleep(5)
      try:
        salome.salome_init(theStudyId=1) #get the number Id 1, or create it
        self.salomeStudy = salome.myStudy
        logger.info("GetSalomeStudy: got study: %s Id: %s",self.salomeStudy._get_Name(),self.salomeStudy._get_StudyId())
        return
      except:
        self.salomeStudy = None
        pass
    else: #except:
      logger.error("I can't get a salome cession. Launch one please: '.../cassis -k'")

if __name__ == '__main__':
  from salomepy.onceQApplication import OnceQApplication
  app = OnceQApplication([''])
  fen = QMainWindowForLog()
  fen.show()
  app.exec_()

