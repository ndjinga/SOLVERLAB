#!/usr/bin/env python
# -*- coding: utf-8 -*-

# %% LICENSE_SALOME_CEA_BEGIN
# see PACKAGESPY/LICENCE file
# %% LICENSE_END


from PyQt5 import QtGui, QtCore, QtWidgets
import os
import sys
import shutil
import fnmatch
import glob
from time import sleep
from datetime import datetime

import salomepy.iconsUser as IUSR
from xyzpy.guiXyz.treeXmlXyz import TreeXmlXyzMainWidget
from salomepy.qTabMultipleTextCentral import QTabMultipleTextCentral
from matplotlibpy.controllerPlt import ControllerPlt
from pandaspy.pandasOscarMainWidget import PandasMainWidget
from matplotlibpy.matplotlibWindowToolbar import MatplotlibWindowToolbar

import xyzpy.loggingXyz as LOG

logger = LOG.getLogger()

verbose = False
verboseEvent = False


class QMainWindowForOscar(QtWidgets.QMainWindow):
  index = [0]

  def __init__(self, *args, **kwargs):
    super(QMainWindowForOscar, self).__init__(*args)
    self.setObjectName("QMainWindowForOscar"+str(self.index))
    self.index[0] += 1 #unambigous objectName
    self.setWindowTitle(self.objectName())
    self.setWindowModality(QtCore.Qt.NonModal)
    self.prefixShortcut = "Ctrl+"
    self.treeViewXml = None
    self._controllerPlt = ControllerPlt() #without desktop, only used to get widgets and toolbar
    self._pandasOscar = PandasMainWidget(plotView = self._controllerPlt.centralPltView) #used to get widgets and toolbar
    self._windowForLog = None
    self.__addCentral()
    self.__addDocks()
    self.__createActions()
    self.__addToolBars()
    self.args = args
    self.kwargs = kwargs
    #print "args, kwargs", args, kwargs 
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
    #self.launchOnArgs.connect(self.theLaunchOnArgs)
    #self.treeViewXml.quickEditFiles.connect(self.centralWidget().quickEditFiles)
    #self.centralWidget().logCastemWidget.finishPopen.connect(self.treeView.refreshModel)
    self.initializeWorkdirDone = False
    self.initializeWorkdir()
    #self.controller.connectTreeViewAndTableEdit(self.treeView, self.centralWidget().csvEditWidget)
    self._controller = None
    self.setTabPosition(QtCore.Qt.LeftDockWidgetArea, QtWidgets.QTabWidget.North)
    for dock in self.docks[1:]:
      self.tabifyDockWidget(self.docks[0], dock)
    if "withDocks" in self.kwargs: 
      withDocks = bool(self.kwargs["withDocks"])
    else:
      withDocks = True
    if withDocks == True:
      for dock in self.docks: dock.show()
    else:
      for dock in self.docks: dock.hide()

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
    central = self._pandasOscar.widget #QTabMultipleTextCentral()
    self.setCentralWidget(central)
    #self.centralWidget().resize(self.centralWidget().size())
    self.centralWidget().show()

  def __addDocks(self):
    self.docks = []
    dock = QtWidgets.QDockWidget("XmlObject", self)
    self.treeViewXml=TreeXmlXyzMainWidget()
    logger.debug("create treeViewXml")
    dock.setWidget(self.treeViewXml.widget)
    dock.setAllowedAreas(QtCore.Qt.LeftDockWidgetArea)
    self.docks.append(dock)
    for dock in self._controllerPlt.docks:
      logger.warning("add dock of _controllerPlt")
      dock.setAllowedAreas(QtCore.Qt.LeftDockWidgetArea)
      self.docks.append(dock)
    for dock in self.docks:
      self.addDockWidget(QtCore.Qt.LeftDockWidgetArea, dock)

  def __addToolBars(self):
    self.toolBars = []
    for tb in self.treeViewXml.toolBars:
      self.toolBars.append(tb)
      self.addToolBar(tb)
    for tb in self._pandasOscar.toolBars:
      self.toolBars.append(tb)
      self.addToolBar(tb)
    for tb in self._controllerPlt.toolBars:
      self.toolBars.append(tb)
      self.addToolBar(tb)

    tb = self.addToolBar("EditLog")
    for action in self.actions:
      tb.addAction(action)
    self.toolBars.append(tb)
    
  def __createActions(self):
    self.actions = []
    #self.actions.append(self.__createAction( "Dgibi", "D", "Show file dgibi", self.showDgibiAction))
    #self.actions.append(self.__createAction( "Castem", "C", "Launch Castem calculus", self.castemLaunchAction, "castem"))
    #self.actions.append(self.__createAction( "Paraview", "S", "Launch Paraview pattern", self.showResultsAction, "paraview"))
    self.actions.append(self.__createAction( "AllTestLauncherOscar", "T", "Launch all python tests", self.testAllAction, "tests"))
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
      logger.warning("create inexisting $WORKDIR4LOG '%s'" % workDir)
      os.makedirs(workDir)
    logger.debug("QMainWindowForOscar logs files in $WORKDIR4LOG '%s'" % workDir)
    #self.treeView.refresh.emit()
    
  def med2sauv(self, file_in, file_out=None):
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

  def launchCmdIntoPopen(self, cmd, clearBefore=True):
    try:
      self._windowForLog.show()
    except:
      if True: #try:
        from salomepy.qMainWindowForLog import QMainWindowForLog
        self._windowForLog = QMainWindowForLog(withDocks=False)
        self._windowForLog.setElementaryLog()
        self._windowForLog.setWindowTitle("Launch Command")
        self._windowForLog.show()
      else: #except:
        logger.error("No window for log, not executed command '%s'"%cmd)
        return
    self._windowForLog.launchCmdIntoPopen("echo coucou", clearBefore=clearBefore)
  
  def testAllAction(self):
    logger.info("testAllAction")
    """ 
    import distenepy.utilsDisteneHelp as UDIS
    data = UDIS.hybridHelp()
    self.treeView.setFromXml(data.toXml())
    self.centralWidget().xmlWidget.insertText(data.toStrXml())
    """
    #unittests, for big example, AllTestLauncher.sh have to be in $PATH
    cmd = "cd $PACKAGESPY_ROOT_DIR; pwd ; AllTestLauncher.sh"
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
  fen = QMainWindowForOscar()
  fen.show()
  app.exec_()

