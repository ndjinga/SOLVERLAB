#!/usr/bin/env python
# -*- coding: utf-8 -*-

# %% LICENSE_SALOME_CEA_BEGIN
# see PACKAGESPY/LICENCE file
# %% LICENSE_END


import os
from PyQt5 import QtGui, QtCore, QtWidgets
import xml.etree.ElementTree as ET
from datetime import datetime

from matplotlibpy.qMainWindowForPlt import QMainWindowForPlt
from matplotlibpy.treeXmlPlt import TreeXmlPlt

import xyzpy.loggingXyz as LOG
from xyzpy.controllerXyz import ControllerXyz
from xyzpy.baseXyz import _XyzConstrainBase, BaseFreeXyz #, ListOfBaseXyz
import xyzpy.baseXyz as BXYZ
import xyzpy.actionsFactoryXyz as ACFX
import xyzpy.intFloatListXyz as IFLX #append factory classes
#from matplotlibpy.MatplotlibPlt import MatplotlibPlt #append factory classes

import xyzpy.utilsXyz as UXYZ
from salomepy.browserPythonClass import BrowserPythonClass


logger = LOG.getLogger()
verbose = False
verboseEvent = False


########################################################################################
class FigurePlt(_XyzConstrainBase):
  """
  as matplotlib Figure, for treeViewXyz, object is hidden in self._figure
  define attribute name and ident for treeView string reference
  """
  _attributesList = [ #list, not a dict because sequential order list is used in files Xyz
    ("name", "StrXyz"),
    ("ident", "StrXyz"),
  ]
  index = [0]

  def __init__(self):
    super(FigurePlt, self).__init__()
    self._figure = None
    self.setIsCast(True)
    self.index[0] += 1 #unambigous objectName
    self._no = self.index[0]
    name = "%s_%i" % (self.getNameObject(), self.index[0]) #could be changed
    ident = "%s_%i" % (self.getDateTimeNow(), self.index[0]) #not to be changed
    self.name = name #FigurePlt_xx
    self.ident = ident #yymmdd_hhmmss_xx
    nameObject = "%s_%s" %(self.getNameObject(), self.ident)
    self.setNameObject(name) #FigurePlt_yymmdd_hhmmss_xx
    self.setIsSet(True)

  def setFigure(self, figure):
    """set 'as it' as hidden attribute: not for toXml method etc..."""
    self._figure = figure

  def getDateTimeNow(self):
    """get an str(), based on datetime.now()"""
    now = datetime.now()
    return now.strftime("%y%m%d_%H%M%S")

########################################################################################
class ModelPlt(BaseFreeXyz):
  """
  firstly a simple set of matplotlib Figure,
  secondly in future other sets of artist ...or else...
  """
  def __init__(self):
    super(ModelPlt, self).__init__()
    self._addOrder = {} #dict, not list, because delete elements allowed
    self._addIndex = 0

  def addFigure(self, aFigure):
    nameFig = str(aFigure.name) #ident is sure, but it is not user_friendly
    if verbose: print("ModelPlt.addFigure %s" % nameFig)
    setattr(self, nameFig, aFigure)
    #print "ModelPlt.addFigure %s\n" % aFigure.name, self.toStrXml()
    self._addOrder[self._addIndex] = nameFig
    self._addIndex += 1

  def deleteFigure(self, name):
    if verbose: print("ModelPlt.deleteFigure %s" % name)
    try:
      delattr(self, name)
    except:
      print("ModelPlt.deleteFigure problem %s\n" % aFigure.name, self.toStrXml())
      return (False, "ModelPlt.deleteFigure problem '%s' not in\n%s" % (name, str(self._addOrder)))
    indexes = [index for index, value in list(self._addOrder.items()) if value==name]
    for index in indexes: del self._addOrder[index]
    #print "delete self._addOrder",self._addOrder, "\n", self
    #for i in self: print "figure",i
    return (True, "")

  def getFigureByName(self, name):
    try:
      res = getattr(self, name)
      return res._figure
    except:
      return None

  def getFigureByIdent(self, ident):
    for index in sorted(self._addOrder.keys()):
      try:
        name = self._addOrder[index]
        res = getattr(self, name)
        if res.ident == ident:
          return res._figure
      except:
        pass
    return None

  def getFigureByNameOrIdent(self, nameOrIdent):
    res = self.getFigureByName(nameOrIdent)
    if res == None: res = self.getFigureByIdent(nameOrIdent)
    return res

  def __iter__(self):
    """iter is on order through _addOrder and so for xml outputs"""
    #print "ModelPlt.__iter__", self.__dict__.keys()
    done = []
    for no in sorted(self._addOrder.keys()):
      name = self._addOrder[no]
      if name not in done and name in list(self.__dict__.keys()):
        done.append(name)
        yield getattr(self, name)
      else:
        print("ModelPlt.__iter__ with __dict__.keys()", list(self.__dict__.keys()))
        raise Exception(self._className + ": attribute '" + name + "' is undefined")


########################################################################################
class ControllerPlt(ControllerXyz):
  """
  class for manage request and action to/from views and model of Matplotlib 1.4
  as MVC pattern
  request data through QtCore.pyqtSignal signals events, are strings Xml
  one controller for one model of some matplotlib Figures and some views
  """

  #http://pyqt.sourceforge.net/Docs/PyQt5/new_style_signals_slots.html
  #http://qt-project.org/doc/qt-5.1/qtwidgets/qaction.html ... Detailed Description
  #LaunchMatplotlibSignal = QtCore.pyqtSignal(object)
  LaunchAllTestsSignal = QtCore.pyqtSignal(object)
  #LaunchMatplotlibCodeTestsSignal = QtCore.pyqtSignal(object)
  #LoadMatplotlibModelXmlSignal = QtCore.pyqtSignal(object)
  #LoadMatplotlibModelPltSignal = QtCore.pyqtSignal(object)
  #SaveMatplotlibModelXmlSignal = QtCore.pyqtSignal(object)
  #SaveMatplotlibModelPltSignal = QtCore.pyqtSignal(object)
  RefreshMatplotlibModelSignal = QtCore.pyqtSignal(object)
  ClearModelSignal = QtCore.pyqtSignal(object)

  def __init__(self, *args, **kwargs):
    super(ControllerPlt, self).__init__(*args, **kwargs)
    self.initializeDone = False
    self.workDir = None     #will not change
    self.currentDir = None  #could change, subdirectory of workDir (a priori)
    self.centralPltView = None
    self.MatplotlibNameFileSaveXml = None
    self.MatplotlibNameDirSavePlt = None
    self.actions = []
    self.treeViews = []
    self.docks = []
    self.toolBars = []
    self.__initialize()
    self._browsers = []
    self._model = None
    self._pltViews = []

  def addNewFigure(self, figure):
    """get new matplotlib figure in treeView browser"""
    aFigure = FigurePlt()
    aFigure.setFigure(figure)
    if verbose: print("ControllerPlt.addNewFigure", aFigure.name)
    """
    if len(self._browsers) == 0: #do once for debug
      #res = BrowserPythonClass(figure, levelMax=1)
      #self._browsers.append(res.display())
      res = BrowserPythonClass(aFigure, levelMax=2)
      self._browsers.append(res.display())
    """
    if self._model == None: self._model = ModelPlt()
    self._model.addFigure(aFigure)
    self.RefreshMatplotlibModelSignal.emit(None)

  def deletePlt(self, namePlt):
    res = self._model.deleteFigure(namePlt)
    self.RefreshMatplotlibModelSignal.emit(None)
    return res

  def displayPlt(self, namePlt):
    aFigure = self._model.getFigureByNameOrIdent(namePlt)
    if aFigure == None: return (False, "no plot for %s" % namePlt)
    aPltView = self.centralPltView
    aPltView.setFigure(aFigure) #new standalone window
    aPltView.setWindowTitle(namePlt)
    return (True,  "ok")

  def displayStandalonePlt(self, namePlt):
    """display in new standalone window"""
    aFigure = self._model.getFigureByNameOrIdent(namePlt)
    if aFigure == None: return (False, "no plot for %s" % namePlt)
    aPltView = QMainWindowForPlt()
    aPltView.setFigure(aFigure) #new standalone window
    aPltView.setWindowTitle(namePlt)
    self._pltViews.append(aPltView)
    aPltView.show()
    return (True,  "ok")


  def __initialize(self):
    if self.initializeDone == True:
      logger.error("initialize only once for controller %s" % (self.objectName()))
      return
    self.__initializeWorkdir()
    self.__createActions()
    self.__addToolBars()
    self.__addCentral()
    self.__addDocks()
    if self._desktop: #salome?
      if verbose: print("ControllerPlt: there is desktop set dock and central widget", self._desktop)
      for dock in self.docks:
        self._desktop.addDockWidget(QtCore.Qt.LeftDockWidgetArea, dock)
        dock.show()
      self._desktop.setCentralWidget(self.centralPltView)
      self._desktop.centralWidget().show()
      for tb in self.toolBars:
        self._desktop.addToolBar(tb)

  def __addCentral(self):
    self.centralPltView = QMainWindowForPlt()
    self.centralPltView.setController(self)
    if self._desktop: #salome? no need dock of QMainWindowForPlt
      for dock in self.centralPltView.docks: dock.hide() #no need
      for tb in self.centralPltView.toolBars: tb.hide()
    logger.debug("create central view for log Matplotlib")
    self.setView(self.centralPltView)
    """
    self.setCentralWidget(central)
    self.centralWidget().resize(self.centralWidget().size())
    self.centralWidget().show
    """

  def __addDocks(self):
    #self.docks = []
    dock = QtWidgets.QDockWidget("MatplotlibObjects")
    treeView = TreeXmlPlt()
    self.treeViews.append(treeView)
    logger.debug("create treeView for Matplotlib")
    dock.setWidget(treeView)
    dock.setAllowedAreas(QtCore.Qt.LeftDockWidgetArea)
    self.docks.append(dock)
    self.setView(treeView) #contains treeViews
    """
    for dock in self.docks:
      self.addDockWidget(QtCore.Qt.LeftDockWidgetArea, dock)
      dock.show()
    """

  def __addToolBars(self):
    #self.toolBars = []
    tb = QtWidgets.QToolBar("EditPlt")      #self.addToolBar("Edit")
    for action in self.actions:
      tb.addAction(action)
    act = ACFX.getCommonActionByName("GeneralHelp")
    if act != None: tb.addAction(act)
    self.toolBars.append(tb)

  def __createActions(self):
    """create actions for self widget AND other widgets through ACFX.addInCommonActions"""
    logger.debug("create actions %s" % (self.objectName()))
    #self.actions = []
    """action = ACFX.QActionXyz(name="LoadMatplotlibModelXml", text="Load Matplotlib data xml")
    ok = action.setAction( slot=self.LoadMatplotlibModelXmlAction, signal=self.LoadMatplotlibModelXmlSignal,
                           shortcut=None, tooltip=u"Load Matplotlib data from file xml", icon="openxml" )
    if ok:
      ACFX.addInCommonActions(action)
      self.actions.append(action)

    action = ACFX.QActionXyz(name="LoadMatplotlibModelPlt", text="Load Matplotlib data cnf")
    ok = action.setAction( slot=self.LoadMatplotlibModelPltAction, signal=self.LoadMatplotlibModelPltSignal,
                           shortcut=None, tooltip=u"Load Matplotlib data from files cnf in directory", icon="opencnf" )
    if ok:
      ACFX.addInCommonActions(action)
      self.actions.append(action)

    action = ACFX.QActionXyz(name="SaveMatplotlibModelXml", text="Save Matplotlib data Xml format")
    ok = action.setAction( slot=self.SaveMatplotlibModelXmlAction, signal=self.SaveMatplotlibModelXmlSignal,
                           shortcut=None, tooltip=u"Save Matplotlib data to file xml", icon="savexml" )
    if ok:
      ACFX.addInCommonActions(action)
      self.actions.append(action)

    action = ACFX.QActionXyz(name="SaveMatplotlibModelPlt", text="Save Matplotlib Data cnf format")
    ok = action.setAction( slot=self.SaveMatplotlibModelPltAction, signal=self.SaveMatplotlibModelPltSignal,
                           shortcut=None, tooltip=u"Save Matplotlib data to files cnf in directory", icon="savecnf" )
    if ok:
      ACFX.addInCommonActions(action)
      self.actions.append(action)

    action = ACFX.QActionXyz(name="LaunchMatplotlib", text="Launch Matplotlib calculus")
    ok = action.setAction( slot=self.LaunchMatplotlibAction, signal=self.LaunchMatplotlibSignal,
                           shortcut=None, tooltip=u"Launch Matplotlib calculus", icon="run" )
    if ok:
      ACFX.addInCommonActions(action)
      self.actions.append(action)
    """

    """action = ACFX.QActionXyz(name="LaunchMatplotlibCodeTests", text="Launch MatplotlibCode tests")
    ok = action.setAction( slot=self.LaunchMatplotlibCodeTestsAction, signal=self.LaunchMatplotlibCodeTestsSignal,
                           shortcut=None, tooltip=u"Launch all MatplotlibCode tests", icon="testMatplotlibCode" )
    if ok:
      ACFX.addInCommonActions(action)
      self.actions.append(action)
    """

    action = ACFX.QActionXyz(name="RefreshMatplotlibModel", text="Refresh views")
    ok = action.setAction( slot=self.RefreshMatplotlibModelAction, signal=self.RefreshMatplotlibModelSignal,
                           shortcut=None, tooltip="Refresh Matplotlib model views", icon="refresh" )
    if ok:
      ACFX.addInCommonActions(action) #may be not in common, only for controller?
      self.actions.append(action)

    action = ACFX.QActionXyz(name="ClearModel", text="Clear Model")
    ok = action.setAction( slot=self.ClearModelAction, signal=self.ClearModelSignal,
                           shortcut=None, tooltip="Clear Matplotlib model", icon="clearModel" )
    if ok:
      ACFX.addInCommonActions(action) #may be not in common, only for controller?
      self.actions.append(action)

    action = ACFX.QActionXyz(name="LaunchAllMatplotlibTests", text="Launch tests")
    ok = action.setAction( slot=self.LaunchAllTestsAction, signal=self.LaunchAllTestsSignal,
                           shortcut=None, tooltip="Launch all Matplotlib tests", icon="test" )
    if ok:
      ACFX.addInCommonActions(action)
      self.actions.append(action)


  def __initializeWorkdir(self):
    """
    initialize user working directory $MATPLOTLIB_WORKDIR if not existing
    """
    nameVar = "MATPLOTLIB_WORKDIR"
    workDir = os.getenv(nameVar)
    if workDir == None:
      homeDir = os.getenv("DEFAULT_WORKDIR", "${HOME}/WORKDIR")
      workDir = os.path.join(homeDir, "MATPLOTLIB")
    workDir = os.path.realpath(os.path.expandvars(workDir))
    os.environ[nameVar] = workDir
    if not os.path.exists(workDir):
      os.makedirs(workDir)
    self.workDir = workDir
    self.currentDir = workDir

  ############################################################################
  """
  slots:
  all name beginning with uppercase are slot without argument
  they are also postfixed ... xxxxAction
  return status is only to significate method is ok, but nothing about action
  """
  ############################################################################
  def ShowAllViewsAction(self):
    if verboseEvent: print("ShowAllViewsAction")
    for view in self.getViews():
      view.show()
    return True

  def ClearModelAction(self):
    if verboseEvent: print("ClearModelAction")
    self.clearModel()
    return True

  def RefreshMatplotlibModelAction(self):
    """refresh Views"""
    if verboseEvent: print("RefreshMatplotlibModelAction")
    #aDataXml =  self._model.toXml() #as a copy
    for treeView in self.treeViews:
      #treeView.update()
      treeView.refreshModel()
    return True

  def LaunchAllTestsAction(self):
    self.centralPltView.allTests()
    #cmd = cmd = "cd $PACKAGESPY_ROOT_DIR; pwd ; AllTestLauncher.sh"
    #if verboseEvent: print "LaunchAllTestsAction '%s'" % cmd
    #self.centralPltView.launchCmdIntoPopen(cmd)
    return True

  def LaunchMatplotlibCodeTestsAction(self):
    cmd = "AllMatplotlibCodeTestLauncher.sh"
    if verboseEvent: print("LaunchMatplotlibCodeTestsAction '%s'" % cmd)
    self.centralPltView.launchCmdIntoPopen(cmd)
    return True

  def LaunchMatplotlibAction(self):
    if verboseEvent: print("LaunchMatplotlibAction")
    #cmd = "AllTestLauncher.py 2>&1" #unittests strerr to stdout
    res = self.SaveMatplotlibModelPltAction()
    if res == False: return
    cmd = "cd %s; MatplotlibCode.sh" % self.MatplotlibNameDirSavePlt
    self.centralPltView.launchCmdIntoPopen(cmd)
    return True

  def SaveMatplotlibModelXmlAction(self):
    if verboseEvent: print("SaveMatplotlibModelXml")
    if self._model == None:
      QtWidgets.QMessageBox.warning(self._desktop, "warning", "no Matplotlib data to save")
      return False
    self.MatplotlibNameFileSaveXml = None #every saveAs
    if self.MatplotlibNameFileSaveXml == None:
      nameFile = QtWidgets.QFileDialog.getSaveFileName(self._desktop, 'Save File', self.currentDir)[0]
    else:
      nameFile = self.MatplotlibNameFileSaveXml
    nameFile = str(nameFile)
    if nameFile == "": return True #cancel
    if ".xml" != nameFile[-4:] : nameFile += ".xml"
    realPath = os.path.realpath(nameFile)
    dirName = os.path.dirname(realPath)
    if not os.path.exists(dirName): os.makedirs(dirName)
    self._model.toFileXml(realPath)
    QtWidgets.QMessageBox.information(self._desktop, "info", "saved file:\n"+realPath)
    self.MatplotlibNameFileSaveXml = realPath #permit resave... but not yet
    return True

  def LoadMatplotlibModelXmlAction(self):
    from widgetpy.salomeQFileDialog import SalomeQFileDialog
    if verboseEvent: print("LoadMatplotlibModelXml")
    aDialog = SalomeQFileDialog()
    nameDir = aDialog.browseDirDialog('Load Microgen data from directory', self.currentDir, ["MATPLOTLIB_WORKDIR", "EXAMPLESMATPLOTLIB"] )
    if nameFile == "": return True #cancel
    realPath = os.path.realpath(nameFile)
    if not os.path.isfile(realPath):
      QtWidgets.QMessageBox.warning(self._desktop, "warning",
          "Load Matplotlib data: not a file \n'%s'" % realPath)
      return False
    try:
      aData = UXYZ.fromFileXml(realPath)
      if aData.__class__.__name__ != "MatplotlibPlt":
        QtWidgets.QMessageBox.warning(self._desktop, "warning", \
            "Load Matplotlib data: '%s' is not a MatplotlibPlt instance from \n'%s'" % \
            (aData.__class__.__name__, realPath) )
        return False
      self.clearModel()
      self.setModel(aData)
      self.RefreshMatplotlibModelSignal.emit(None)
      return True
    except Exception as e:
      QtWidgets.QMessageBox.warning(self._desktop, "warning",
          "Load Matplotlib data: problem loading\n'%s'" % e)
      return False

if __name__ == '__main__':
    from salomepy.onceQApplication import OnceQApplication
    app = OnceQApplication()
    desktop = QtWidgets.QMainWindow()
    desktop.resize(1000, 600)
    ctrl = ControllerPlt(desktop = desktop)
    desktop.setWindowTitle(ctrl.objectName())
    desktop.show()
    app.exec_()
