#!/usr/bin/env python
# -*- coding: utf-8 -*-

# %% LICENSE_SALOME_CEA_BEGIN
# Copyright (C) 2008-2023  CEA
#
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 2.1 of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with this library; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
#
# See http://www.salome-platform.org or email : webmaster.salome@opencascade.com
# %% LICENSE_END


import os
import platform
import pprint
import time
import shutil
import pprint as PP
import subprocess as SP

from PyQt5 import QtGui, QtCore, QtWidgets as QTW

import xyzpy.loggingXyz as LOG
from xyzpy.controllerXyz import ControllerXyz
import xyzpy.actionsFactoryXyz as ACFX
import debogpy.debug as DBG

# append factory classes
import xyzpy.intFloatListXyz as IFLX  # append factory classes

import solverlabpy.modelSvl
from solverlabpy.modelSvl import ModelSvl

from salomepy.qTabMultipleTextCentral import QTabMultipleTextCentral
from salomepy.threadWorkerForWidgetEdit import ThreadWorkerForWidgetEdit
import sys
# from qtRootCanvasWidget import QtRootCanvasWidget

import xyzpy.utilsXyz as UXYZ
import salomepy.utilsWorkdir as UTW
import solverlabpy.solverlabFilePatterns as UFPA
import solverlabpy.solverlabSettings as USET

import solverlabpy.LaunchDiffusionEq as LDEQ
logger = LOG.getLogger()
verbose = True
verboseEvent = True

_MyDir = os.path.split(os.path.realpath(__file__))[0]

try:
  import solverlab
except Exception as e:
  logger.warning("problem import solverlab : %s" % e)


##############################################################################################################################################
def join(*v):
  """
  as os.path.join but set antislash as slash,
  even for windows, keep windows 'c:'
  """
  # https://stackoverflow.com/questions/12086224/why-not-os-path-join-use-os-path-sep-or-os-sep
  res = os.path.join(*v).replace("\\", "/")
  return res


########################################################################################
class QTabMultipleTextCentralSvl(QTabMultipleTextCentral):
  # -1 as useless as not created
  # user use inheritage for modify TAB_... order or appearence tabs
  TAB_LOG_CMD = 0
  TAB_OTHERTEXTEDIT = 2
  TAB_FILESYSTEM = 1

  # additional inherited
  # TAB_ROOTCANVAS = 3

  def _init_tabs_inherited(self):
    return
    """
    # case usage of ROOT, and URANIE, for future
    if self.TAB_ROOTCANVAS != -1:
      wid = QtRootCanvasWidget()
      wid.saveFileExt = ".solverlab"
      wid.tabName = "Canvas"  # as foreign clear short label name
      # wid.attName = "ROOTCanvasWidget"  # as immmutable name
      # self.activeTabs[self.TAB_ROOTCANVAS] = wid
    """


########################################################################################
class ControllerSvl(ControllerXyz):
  """
  class for manage request and action to/from views and model of solverlabGUI
  as MVC pattern
  one controller for one model and some views
  """

  # http://pyqt.sourceforge.net/Docs/PyQt4/new_style_signals_slots.html
  # http://qt-project.org/doc/qt-5.1/qtwidgets/qaction.html ... Detailed Description

  LaunchSolverlabSignal = QtCore.pyqtSignal(object)
  LaunchAllTestsSignal = QtCore.pyqtSignal(object)
  LaunchSolverlabCodeTestsSignal = QtCore.pyqtSignal(object)
  LoadSolverlabDefaultModelSignal = QtCore.pyqtSignal(object)
  LoadSolverlabModelXmlSignal = QtCore.pyqtSignal(object)
  SaveSolverlabModelXmlSignal = QtCore.pyqtSignal(object)
  RefreshSolverlabModelSignal = QtCore.pyqtSignal(object)
  LoadDataFileInModelSignal = QtCore.pyqtSignal(object)
  UpdateEtudeSignal = QtCore.pyqtSignal(object)
  UpdateEtudeDataSignal = QtCore.pyqtSignal(object)
  UpdateEtudeRootlogonSignal = QtCore.pyqtSignal(object)
  CreateCpackEtudeSignal = QtCore.pyqtSignal(object)
  CreateDocSignal = QtCore.pyqtSignal(object)
  GitCommitEtudeSignal = QtCore.pyqtSignal(object)
  ClearModelSignal = QtCore.pyqtSignal(object)
  ExecPythonCodeSignal = QtCore.pyqtSignal(object)
  ExecRootCodeSignal = QtCore.pyqtSignal(object)
  ExecRootProcessFileSignal = QtCore.pyqtSignal(object)
  ExecSearchMethodInSolverlabSignal = QtCore.pyqtSignal(object)
  ExecPrintROOTContextSignal = QtCore.pyqtSignal(object)
  ExecPostTreatmentsSignal = QtCore.pyqtSignal(object)
  UserExpandSignal = QtCore.pyqtSignal(list)
  UserModeSignal = QtCore.pyqtSignal(str)
  SolverlabGuiHelpSignal = QtCore.pyqtSignal(object)
  SolverlabCodeHelpSignal = QtCore.pyqtSignal(object)
  SolverlabExampleSignal = QtCore.pyqtSignal(object)

  def __init__(self, *args, **kwargs):
    super(ControllerSvl, self).__init__(*args, **kwargs)
    self.initializeDone = False
    self.workDir = None  # will not change
    self.currentDir = None  # could change, subdirectory of workDir (a priori)
    self.centralLogView = None
    self.solverlabNameFileSaveXml = None
    self.solverlabNameDirSaveSvl = None
    self.actions = []
    self.treeViews = []
    self.docks = []
    self.toolBars = []
    self.__initialize()
    self._solverlabXmlDefaultName = "solverlabGui.xml"
    # self.setIpcController(ControllerIpcSvl())

    # connect
    self.UpdateEtudeSignal.connect(self._updateEtude)
    self.UpdateEtudeDataSignal.connect(self._updateDataEtude)
    #self.UpdateEtudeRootlogonSignal.connect(self._updateRootlogonEtude)
    #self.GitCommitEtudeSignal.connect(self._gitCommitEtude)
    #self.CreateCpackEtudeSignal.connect(self._createCpackEtude)
    #self.CreateDocSignal.connect(self._createDocEtude)
    #self.ExecPythonCodeSignal.connect(self._ExecPythonCode)
    #self.ExecRootCodeSignal.connect(self._ExecRootCode)
    #self.ExecRootProcessFileSignal.connect(self._ExecRootProcessFile)
    #self.ExecSearchMethodInSolverlabSignal.connect(self._ExecSearchMethodInSolverlab)
    # self.ExecPrintROOTContextSignal.connect(self._ExecPrintROOTContext)
    self.ExecPostTreatmentsSignal.connect(self._ExecPostTreatments)
    self.UserExpandSignal.connect(self._UserExpand)
    self.UserModeSignal.connect(self._UserMode)
    self.LoadDataFileInModelSignal.connect(self._loadDataFileInModel)
    #self.modelChangeSignal.connect(self.SetSolverlabGUIDIR)

    self.isController = True

    self.SOLVERLABGUI_ROOT_DIR = os.getenv("SOLVERLABGUI_ROOT_DIR")
    self.SOLVERLABGUI_WORKDIR = os.getenv("SOLVERLABGUI_WORKDIR")

  def __initialize(self):
    if self.initializeDone == True:
      logger.error("initialize only once for controller %s" % (self.objectName()))
      return
    self.__initializeWorkdir()
    self.__createActions()
    self.__addToolBars()
    self.__addCentral()
    self.__addDocks()
    if self._desktop:  # salome?
      logger.debug("ControllerSvl: desktop have set dock and central widget")
      for dock in self.docks:
        self._desktop.addDockWidget(QtCore.Qt.LeftDockWidgetArea, dock)
        dock.show()
      self._desktop.setCentralWidget(self.centralLogView)
      self._desktop.centralWidget().show()
      for tb in self.toolBars:
        self._desktop.addToolBar(tb)

  def __addCentral(self):
    import salomepy.qMainWindowForLog as QMFL

    self.centralLogView = QMFL.QMainWindowForLog(centralWidget=QTabMultipleTextCentralSvl())
    if self._desktop:  # salome? no need dock of QMainWindowForLog
      for dock in self.centralLogView.docks: dock.hide()  # no need
      for tb in self.centralLogView.toolBars: tb.hide()

    logger.info("create central view for log solverlab")
    self.setView(self.centralLogView)

    tabs = self.centralLogView.getTabs()
    filetab = self.centralLogView.getTabByName("Explore Dir")
    rootPath = os.path.join("${SOLVERLABGUI_WORKDIR}")
    filetab.setDirRootPath(rootPath, filters=[])
    """
    self.setCentralWidget(central)
    self.centralWidget().resize(self.centralWidget().size())
    self.centralWidget().show
    """

  def __addDocks(self):
    import solverlabpy.treeViewSvl as TVSVL
    # self.docks = []
    dock = QTW.QDockWidget("SolverlabObjects")
    treeView = TVSVL.TreeViewSvl()
    self.treeViews.append(treeView)
    logger.info("create treeView for solverlab")
    dock.setWidget(treeView)
    dock.setAllowedAreas(QtCore.Qt.LeftDockWidgetArea)
    self.docks.append(dock)
    self.setView(treeView)  # contains treeViews
    """
    for dock in self.docks:
      self.addDockWidget(QtCore.Qt.LeftDockWidgetArea, dock)
      dock.show()
    """

  def __addToolBars(self):
    # self.toolBars = []
    tb = QTW.QToolBar("EditSvl")  # self.addToolBar("Edit")
    for action in self.actions:
      tb.addAction(action)
    # act = ACFX.getCommonActionByName("GeneralHelp")
    # if act != None: tb.addAction(act)
    self.toolBars.append(tb)

  def __createActions(self):
    """create actions for self widget AND other widgets through ACFX.addInCommonActions"""
    logger.debug("create actions %s" % (self.objectName()))
    # self.actions = []
    action = ACFX.QActionXyz(name="LoadSolverlabDefaultModel", text="Load Default Solverlab data")
    ok = action.setAction(slot=self.LoadSolverlabDefaultModelAction, signal=self.LoadSolverlabDefaultModelSignal,
                          shortcut=None, tooltip=u"New Solverlab data", icon="solverlabpy.resources.solverlabgui")
    if ok:
      ACFX.addInCommonActions(action)
      self.actions.append(action)

    action = ACFX.QActionXyz(name="LoadSolverlabModelXml", text="Load Solverlab data xml")
    ok = action.setAction(slot=self.LoadSolverlabModelXmlAction, signal=self.LoadSolverlabModelXmlSignal,
                          shortcut=None, tooltip=u"Load Solverlab data from file xml", icon="openxml")
    if ok:
      ACFX.addInCommonActions(action)
      self.actions.append(action)

    action = ACFX.QActionXyz(name="SaveSolverlabModelXml", text="Save Solverlab data XmlXyz format")
    ok = action.setAction(slot=self.SaveSolverlabModelXmlAction, signal=self.SaveSolverlabModelXmlSignal,
                          shortcut=None, tooltip=u"Save Solverlab data to file xml", icon="savexml")
    if ok:
      ACFX.addInCommonActions(action)
      self.actions.append(action)

    action = ACFX.QActionXyz(name="LaunchSolverlab", text="Launch Solverlab calculus")
    ok = action.setAction(slot=self.LaunchSolverlabAction, signal=self.LaunchSolverlabSignal,
                          shortcut=None, tooltip=u"Launch Solverlab calculus", icon="run")
    if ok:
      ACFX.addInCommonActions(action)
      self.actions.append(action)

    if False: #DBG.isDeveloper():
      action = ACFX.QActionXyz(name="LaunchAllTests", text="Launch tests")
      ok = action.setAction(slot=self.LaunchAllTestsAction, signal=self.LaunchAllTestsSignal,
                            shortcut=None, tooltip=u"Launch all tests", icon="test")
      if ok:
        ACFX.addInCommonActions(action)
        self.actions.append(action)

      action = ACFX.QActionXyz(name="LaunchSolverlabCodeTests", text="Launch SolverlabCode tests")
      ok = action.setAction(slot=self.LaunchSolverlabCodeTestsAction, signal=self.LaunchSolverlabCodeTestsSignal,
                            shortcut=None, tooltip=u"Launch all SolverlabCode tests", icon="testSolverlabCode")
      if False: # ok: 2020 not done yey
        ACFX.addInCommonActions(action)
        self.actions.append(action)

    action = ACFX.QActionXyz(name="RefreshSolverlabModel", text="Refresh views")
    ok = action.setAction(slot=self.RefreshSolverlabModelAction, signal=self.RefreshSolverlabModelSignal,
                          shortcut=None, tooltip=u"Refresh SolverlabObjects tree view", icon="refresh")
    if ok:
      ACFX.addInCommonActions(action)  # may be not in common, only for controller?
      self.actions.append(action)

    action = ACFX.QActionXyz(name="ClearModel", text="Clear Model")
    ok = action.setAction(slot=self.ClearModelAction, signal=self.ClearModelSignal,
                          shortcut=None, tooltip=u"Clear Solverlab model", icon="clearModel")
    if ok:
      ACFX.addInCommonActions(action)  # may be not in common, only for controller?
      self.actions.append(action)

    action = ACFX.QActionXyz(name="GuiHelp", text="solverlab GUI help")
    ok = action.setAction(slot=self.SolverlabGuiHelpAction, signal=self.SolverlabGuiHelpSignal,
                          shortcut=None, tooltip=u"solverlab GUI help", icon="helpGui")
    if ok:
      ACFX.addInCommonActions(action)  # may be not in common, only for controller?
      self.actions.append(action)

    action = ACFX.QActionXyz(name="SolverlabHelp", text="solverlab CODE help")
    ok = action.setAction(slot=self.SolverlabCodeHelpAction, signal=self.SolverlabCodeHelpSignal,
                          shortcut=None, tooltip=u"solverlab CODE help", icon="helpCode")
    if ok:
      ACFX.addInCommonActions(action)  # may be not in common, only for controller?
      self.actions.append(action)

    action = ACFX.QActionXyz(name="SolverlabExample", text="solverlab Example")
    ok = action.setAction(slot=self.SolverlabExampleAction, signal=self.SolverlabExampleSignal,
                          shortcut=None, tooltip=u"load one of .xml to see an example", icon="test")
    if ok:
      ACFX.addInCommonActions(action)
      self.actions.append(action)

  def SolverlabGuiHelpAction(self):
    nameBrowser = UXYZ.getBrowser()
    tmp = os.path.join(self.SOLVERLABGUI_ROOT_DIR, 'doc', 'build', 'html', 'index.html')
    nameUrlHelp = os.path.expandvars(tmp)
    if os.path.exists(nameUrlHelp):
      cmd = "%s %s &" % (nameBrowser, nameUrlHelp)
      proc = SP.Popen(cmd, shell=True)
    else:
      logger.error("inexisting name Url Help: '%s'" % nameUrlHelp )

  def SolverlabCodeHelpAction(self):
    nameBrowser = UXYZ.getBrowser()
    tmp = os.path.join(self.SOLVERLABGUI_ROOT_DIR, 'share', 'doc', 'CoreFlows.pdf')
    nameUrlHelp = os.path.expandvars(tmp)
    if os.path.exists(nameUrlHelp):
      cmd = "%s %s &" % (nameBrowser, nameUrlHelp)
      proc = SP.Popen(cmd, shell=True)
    else:
      logger.error("inexisting name Url Help: '%s'" % nameUrlHelp )


  def __initializeWorkdir(self):
    """
    initialize user working directory $SOLVERLABGUI_WORKDIR if not existing
    """
    nameVar = "SOLVERLABGUI_WORKDIR"
    workDir = os.getenv(nameVar)
    if workDir == None:
      homeDir = os.getenv("HOME")
      workDir = os.path.join(homeDir, "SOLVERLABGUI_WORKDIR")
    workDir = os.path.realpath(workDir)
    os.environ[nameVar] = workDir
    if not os.path.exists(workDir):
      os.makedirs(workDir)
    self.workDir = workDir
    self.currentDir = workDir

  def getEtudeWorkdirExpanded(self):
    if self._model == None:
      return None
    else:
      return self._model.getEtudeWorkdirExpanded()

  ############################################################################
  """
  ############################################################################
  slots:
  all name beginning with uppercase are slot without argument
  they are also postfixed ... xxxxAction
  return status is only to significate method is ok, but nothing about action
  ############################################################################
  """

  def ShowAllViewsAction(self):
    logger.debug("ShowAllViewsAction")
    for view in self.getViews():
      view.show()
    return True

  def ClearModelAction(self):
    logger.debug("ClearModelAction")
    self.clearModel()
    return True

  def RefreshSolverlabModelAction(self):
    """refresh Views"""
    logger.debug("RefreshSolverlabModelAction")
    # aDataXml =  self._model.toXml() #as a copy
    # self.SetSolverlabGUIDIR() #SolverlabGUIDIR change?
    eq = self._model.getSelectedEquation()
    if eq is not None:
      # eq.boundary_condition.update()
      logger.warning("TODO RefreshSolverlabModelAction eq is not None :\nmay be eq.boundary_condition.update")
    for treeView in self.treeViews:
      # treeView.update()
      treeView.refreshModel()
    return True

  def ExpandAllAction(self):
    logger.debug("ExpandAllAction")
    for treeView in self.treeViews:
      treeView.expandAll()
    return True

  def LaunchAllTestsAction(self):
    # cmd = "AllTestLauncher.py 2>&1" #unittests strerr to stdout
    filesh = os.path.join(self.SOLVERLABGUI_ROOT_DIR, 'bin', 'AllSolverlabGuiTestLauncher.sh')
    if os.path.exists(filesh):
      cmd = "cd " + self.SOLVERLABGUI_ROOT_DIR + "; pwd ; ./bin/AllSolverlabGuiTestLauncher.sh"
    else:
        cmd = "cd " + self.SOLVERLABGUI_ROOT_DIR + "; pwd ; ./bin/AllTestLauncher.sh"
    logger.debug("LaunchAllTestsAction '%s'" % cmd)
    self.centralLogView.launchCmdIntoPopen(cmd)
    return True

  def LaunchSolverlabCodeTestsAction(self):
    # 2020 TODO
    cmd = "AllSolverlabCodeTestLauncher.sh"
    logger.debug("LaunchSolverlabCodeTestsAction '%s'" % cmd)
    self.centralLogView.launchCmdIntoPopen(cmd)
    return True

  def SolverlabExampleAction(self):
    from widgetpy.salomeQFileDialog import SalomeQFileDialog
    aDialog = SalomeQFileDialog(parent=self._desktop)
    aDir = os.path.join(self.SOLVERLABGUI_ROOT_DIR, "example", "model")
    nameFile = aDialog.browseFileDialog('Load Solverlab file xml', aDir, "(*.xml *.XML)", [])
    logger.debug("LoadSolverlabModelXml %s" % nameFile)
    if nameFile == "": return True  # cancel
    realPath = os.path.realpath(nameFile)
    if not os.path.isfile(realPath):
      QTW.QMessageBox.warning(self._desktop, "warning",
                              "Load Solverlab xml data: not a file \n'%s'" % realPath)
      return False
    try:
      aData = UXYZ.fromFileXml(realPath)
      logger.debug("LoadSolverlabModelXml: data loaded from file")
      if aData.__class__.__name__ != "ModelSvl":
        QTW.QMessageBox.warning(self._desktop, \
                                "warning", "Load Solverlab data: '%s' is not a ModelSvl instance from \n'%s'" % (
                                aData.__class__.__name__, realPath))
        return False
      self.clearModel()
      self.setModel(aData)
      # logger.debug("LoadSolverlabModelXml: set model data done")
      self.RefreshSolverlabModelSignal.emit(None)
      self._model.userExpand()
      self.ExpandAllAction()
      return True
    except Exception as e:
      # traceback.print_exc()
      QTW.QMessageBox.warning(self._desktop, "warning", "Load Solverlab xml data: problem loading\n%s" % e)
      return False

  def CreateEtudeWorkdir(self):
    etudeDir = self._model.getEtudeWorkdirExpanded()
    etudeDirBrut = self._model.getEtudeWorkdirBrut()
    if os.path.isfile(etudeDir):
      QTW.QMessageBox.warning(self._desktop, "warning",
            "etude directory existing yet as file, fix it:\n'%s'" % etudeDir)
      return False

    if os.path.isdir(etudeDir):
      # QTW.QMessageBox.warning(self._desktop, "warning", "etude directory existing yet:\n'%s'" % etudeDir)
      pass
    else:
      UTW.makeDir(etudeDir)
    # for subDir in "data macros doc".split():
    for subDir in "output".split():
      aDir = os.path.join(etudeDir, subDir)
      if not os.path.isdir(aDir):
        UTW.makeDir(aDir)
        logger.info("create directory %s" % aDir)
    #QTW.QMessageBox.information(self._desktop, "info", "current etude directory:\n'%s'" % etudeDir)
    #self._createGitignore(etudeDir)
    #self._createREADMEfile(etudeDirBrut)
    #self._createDoxyfile(etudeDirBrut)
    #self._createRootLogon(etudeDirBrut)
    #self._gitInit(etudeDir)
    #self._gitCommit("initial commit", etudeDir)
    return True

  def createFilePattern(self, nameFile, patternName):
    etudeDirExp = self._model.getEtudeWorkdirExpanded()
    if not os.path.isdir(etudeDirExp):
      QTW.QMessageBox.warning(self._desktop, "warning", "etude directory not existing, create it.")
      self.CreateEtudeWorkdir()

    replaces = [("@FILE@", nameFile)]
    existingPatterns = UFPA.getPatternKeys()
    if not patternName in existingPatterns:
      _, ext = os.path.splitext(nameFile)
      patternNameByDefault = "aUserFile%s" % ext
      contents = UFPA.getFilePatterns(patternNameByDefault, replaces=replaces)
    else:
      contents = UFPA.getFilePatterns(patternName, replaces=replaces)
    with open(nameFile, "w") as f:
      f.write(contents)
    return True

  def _insertFileInListOf(self, index, aFile, theListOf):
    """no insert if existing yet"""
    for i in theListOf:
      if str(i) == aFile: return True  # existing yet TODO compare true expandvars name file
    try:
      newItem = theListOf._allowedClasses[0](aFile)  # immutable
    except:
      newItem = theListOf._allowedClasses[0]()  # mutable as data List Of
      newItem.name = aFile
    if index >= 0:
      theListOf.insert(index, newItem)
    else:  # index -1
      theListOf.append(newItem)
    return True

  def _createREADMEfile(self, etudeDirBrut):
    """
    create file README.txt
    """
    name = "README.txt"
    nameFile = os.path.join(etudeDirBrut, name)
    nameFileExp = os.path.expandvars(nameFile)
    res = self.createFilePattern(nameFileExp, name)
    ok = self._insertFileInListOf(0, nameFile, self._model.Analysis.userFileManager)
    self.RefreshSolverlabModelSignal.emit(None)  # changed model
    return True

  def _updateEtude(self, toUpdate="all"):
    """
    override files in etude, with message, and synchronize controller model
    """
    return False
    # TODO obsolete ?
    etudeDirBrut = self._model.getEtudeWorkdirBrut()
    etudeDirExp = self._model.getEtudeWorkdirExpanded()
    if not os.path.isdir(etudeDirExp):
      self.CreateEtudeWorkdir()

    logger.warning("TODO update etude %s create other directories yes or no" % etudeDirBrut)
    return True


    logger.debug("update etude %s" % etudeDirBrut)
    ###### dir functions
    functions = self._model.Analysis.macroManager.functions
    # new functions instance for ultimate replace in model
    arg1 = functions.__class__()
    for i in functions:
      newFile = self._smartCopyFileInEtude(i, "macros")
      arg1.append(newFile)
    cmd = ".Analysis.macroManager.functions = args[1]"
    self.setModelItemValueSignalList.emit([cmd, arg1])

    ###### dir macros
    macros = self._model.Analysis.macroManager.macros
    # new macros instance for ultimate replace in model
    arg1 = macros.__class__()
    for i in macros:
      newFile = self._smartCopyFileInEtude(i, "macros")
      arg1.append(newFile)
    cmd = ".Analysis.macroManager.macros = args[1]"
    self.setModelItemValueSignalList.emit([cmd, arg1])

    ###### dir libraries
    libs = self._model.Analysis.macroManager.libraries
    # new libraries instance for ultimate replace in model
    arg1 = libs.__class__()
    for i in libs:
      newFile = self._smartCopyFileInEtude(i, "macros")
      arg1.append(newFile)
    cmd = ".Analysis.macroManager.libraries = args[1]"
    self.setModelItemValueSignalList.emit([cmd, arg1])
    self.UpdateEtudeRootlogonSignal.emit(None)  # automatic
    return True

  def _loadDataFileInModel(self, aFile):
    return False
    # TODO obsolete ?
    name = aFile
    nameFileExp = os.path.expandvars(name)
    nameFile = os.path.realpath(nameFileExp)
    if not os.path.isfile(nameFile):
      ret = mbox.warning(self._desktop, "error",
                         "inexisting data file:\n%s" % nameFile,
                         mbox.Ok)
      if ret == mbox.Ok: return False  # nothing to do
    theListOf = self._model.Analysis.dataManager
    for datai in theListOf:  # mutable class
      ii = datai.getNameExpanded()  # compare true expandvars name file
      if ii == nameFile:
        self.RefreshSolverlabModelSignal.emit(None)  # no changed model, may be changed file contents?
        self.RefreshSolverlabModelSignal.emit(None)  # changed contents precaution
        return True  # existing yet
    ok = self._insertFileInListOf(-1, name, theListOf)
    # logger.debug("_loadDataFileInModel:\n  %s %s %s" % (aFile, nameFile, PP.pformat(theListOf)))
    self.RefreshSolverlabModelSignal.emit(None)  # changed model
    return True

  def _ExecPostTreatments(self):
    logger.warning("TODO _ExecPostTreatments")
    #from matplotlibpy.matplotlibWindowToolbar import MatplotlibWindowToolbar
    #self.MatplotlibWindow = MatplotlibWindowToolbar()
    #self.MatplotlibWindow.show()

    #from pandaspy.pandasOscarMainWidget import PandasTabWidget, PandasMainWidget
    #from pandaspy.pandasMainWidgetXyz import PandasMainDialogXyz
    import pandaspy.pandasMainWidgetXyz as PDMW
    self.PandasWindow = PDMW.PandasMainDialogXyz("SolverlabGUI Plot Viewer")
    self.PandasWindow.show()
    return True

  def _updateDataEtude(self, toUpdate="all"):
    """
    override files in etude/data, with message, and synchronize controller model
    """
    etudeDirBrut = self._model.getEtudeWorkdirBrut()
    etudeDirExp = self._model.getEtudeWorkdirExpanded()
    if not os.path.isdir(etudeDirExp):
      self.CreateEtudeWorkdir()

    logger.debug("updateDataEtude %s" % etudeDirBrut)

    ###### data files
    dataManager = self._model.Analysis.dataManager
    # new dataManager instance for ultimate replace in model
    arg1 = dataManager.__class__()
    for i in dataManager:
      if toUpdate in ["all", i.name]:
        newFile = self._smartCopyFileDataInEtude(i, "data")
      else:
        newFile = i.__class__();
        newFile.name = str(i.name)
      arg1.append(newFile)

    cmd = ".Analysis.dataManager = args[1]"
    self.setModelItemValueSignalList.emit([cmd, arg1])
    return True

  def _splitFileName(self, aFile):
    """TODO future small tests avoid bug"""
    f = str(aFile)
    aDir, aName = os.path.split(f)  # as ('$HOME', 'aa.txt')
    aDirExp, aNameExp = os.path.split(os.path.realpath(os.path.expandvars(f)))
    return aDir, aName, aDirExp, aNameExp

  def _copyFile(self, originFile, targetFile):
    """copy expanded files and question if override"""
    logger.debug("copy file:\noriginFile %s\ntargetFile %s\n" % (originFile, targetFile))
    if originFile == targetFile: return False  # nothing to do
    if os.path.isfile(targetFile):
      mbox = QTW.QMessageBox
      ret = mbox.warning(self._desktop, "warning",
                         "override file:\n'%s'?" % targetFile,
                         mbox.Yes | mbox.No)
      if ret == mbox.No: return False  # nothing to do
    # copy file...
    shutil.copyfile(originFile, targetFile)
    # message...
    text = "\ncopyFile:\n  originFile %s\n  targetFile %s\n" % (originFile, targetFile)
    self.centralLogView.insertText(text)
    hist = self._model.getHistoryFile()
    hist.appendHistoryCopyOf(originFile, targetFile)
    return True

  def _smartCopyFileInEtude(self, aFile, targetDir):
    """
    have to be smart to copy and override if necessary
    """
    etudeDirBrut = self._model.getEtudeWorkdirBrut()
    etudeDirExp = self._model.getEtudeWorkdirExpanded()

    aDir, aName, aDirExp, aNameExp = self._splitFileName(aFile)

    originFile = os.path.join(aDirExp, aNameExp)
    targetFile = os.path.join(etudeDirExp, targetDir, aNameExp)
    targetFileBrut = os.path.join(etudeDirBrut, targetDir, aNameExp)

    ok = self._copyFile(originFile, targetFile)

    if ok:
      newFile = aFile.__class__(targetFileBrut)
    else:
      newFile = aFile.__class__(str(aFile))  # duplicate reference without copy
    return newFile

  def _smartCopyFileDataInEtude(self, aFile, targetDir):
    """
    have to be smart to copy and override if necessary
    assume dataManager items
    """
    etudeDirBrut = self._model.getEtudeWorkdirBrut()
    etudeDirExp = self._model.getEtudeWorkdirExpanded()

    aDir, aName, aDirExp, aNameExp = self._splitFileName(aFile.name)

    originFile = os.path.join(aDirExp, aNameExp)
    targetFile = os.path.join(etudeDirExp, targetDir, aNameExp)
    targetFileBrut = os.path.join(etudeDirBrut, targetDir, aNameExp)

    ok = self._copyFile(originFile, targetFile)

    if ok:
      newFile = aFile.__class__();
      newFile.name = targetFileBrut
    else:
      newFile = aFile.__class__();
      newFile.name = targetFileBrut  # duplicate reference without copy
    return newFile

  def _launchEtude(self):
    """
    TODO
    """
    #mbox = QTW.QMessageBox
    #ret = mbox.question(self._desktop, "question", " TODO !!! launch background ?", mbox.Yes | mbox.No)
    #logger.warning("ControllerSvl._launchEtude() TODO...")
    self._model.Analysis.toFileSvl()
    return

  def assertDataDirectory(self, etudeDir):
    dataDir = os.path.realpath(os.path.join(etudeDir, "..", "data"))
    if os.path.isdir(dataDir):
      logger.info("use existing corteo data directory\n%s" % dataDir)
    else:
      name = os.path.join(*"${SOLVERLABGUI_ROOT_DIR} solverlabCode data_4bit".split())
      origDir = os.path.realpath(os.path.expandvars(name))
      logger.warning("create user data corteo directory:\n%s ->\n%s" % (origDir, dataDir))
      UTW.copyDir(origDir, dataDir)
  '''
  def runForeground(self):
    stdoutWorker = ThreadWorkerForWidgetEdit(sys.stdout, self, "Black")
    stderrWorker = ThreadWorkerForWidgetEdit(sys.stderr, self, "Red", emit=False)
    # warning 2017 : I never seen finishPopen on stderr don't know why...
    # so set stdout stderr both for garbage collecting on stdoutWorker signal
    log = self.centralLogView.getTabByName("Log Run Code")

    log._Workers[stdoutWorker.objectName()] = (stdoutWorker, stderrWorker)
    if verbose:
      print("%s._Workers %s" % (self.objectName(), PP.pformat(log._Workers)))
    stdoutWorker.start()
    stderrWorker.start()'''

  def LaunchSolverlabAction(self):
    if self._model is None:
      logger.warning("Inexisting solverlab data tree")
      mbox = QTW.QMessageBox
      ret = mbox.question(self._desktop, "question", "Inexisting solverlab data tree,\nCreate one ?", mbox.Yes | mbox.No)
      if ret == mbox.Yes:
        self.LoadSolverlabDefaultModelAction()
      return False

    # launch in foreground
    if not self._model.getBackground():
      LDEQ.launchDEQ(self._model)
      return

    # launch in background
    # create directory ${SolverlabWorkdir}/output/n for the case
    etudeDir = self._model.getEtudeWorkdirExpanded() + "/output"
    i = 1
    while True:
      testdir = etudeDir + "/case%s" % i
      if os.path.isdir(testdir):
        i += 1
      else:
        etudeDir = testdir
        break
    os.makedirs(etudeDir)
    logger.info("Launch solverlab in %s" % etudeDir)

    # ask to save .med in case dir
    mbox = QTW.QMessageBox
    msg = "Do you want to copy the medfile in local ?"
    localcopy = mbox.question(self._desktop, "question", msg, mbox.Yes | mbox.Cancel)
    logger.info("save solverlab model xml")
    if localcopy == mbox.Yes:
      shutil.copy2(self._model.getMedfile(), self.getEtudeWorkdirExpanded())
      name = self._model.GeometryMed.fileMed
      self._model.GeometryMed.fileMed = self.getEtudeWorkdirExpanded() + "/" +os.path.basename(name)
      # save model in case dir with path to local .med
      nameFile = os.path.join(etudeDir, self._solverlabXmlDefaultName)
      self._model.toFileXml(nameFile)
      self._model.GeometryMed.fileMed = name
    else:
      # save model in case dir with path to distant .med
      nameFile = os.path.join(etudeDir, self._solverlabXmlDefaultName)
      self._model.toFileXml(nameFile)

    shutil.copy2(self.SOLVERLABGUI_ROOT_DIR + "/solverlabpy/runSolverlab.py", self.getEtudeWorkdirExpanded())

    if platform.system() == "Windows":
      name = "launchSOLVERLAB.bat"
      nameFile = os.path.join(etudeDir, name)
      strLaunch = self.getStrLaunchBat(etudeDir)
      with open(nameFile, "w") as f:
        f.write(strLaunch)
      # useless .bat windows self.chmodarwx(nameFile)
      cmd = "START /B %s" % nameFile # windows popen accept ONE line
      # cmd = strLaunch
      msg = "Launch solverlab bat script:\n%s" % cmd
    else: #"Linux etc."
      name = "launchSOLVERLAB.sh"
      nameFile = os.path.join(etudeDir, name)
      strLaunch = self.getStrLaunchSh(etudeDir)
      with open(nameFile, "w") as f:
        f.write(strLaunch)
      self.chmodarwx(nameFile)
      cmd = strLaunch # linux popen accept multiples lines
      msg = "Launch solverlab bash script:\n%s" % cmd

    logger.info(msg)
    # mbox = QTW.QMessageBox
    # ret = mbox.question(self._desktop, "question", "Solverlab data files created\nLaunch Solverlab code ?", mbox.Yes | mbox.No)
    # if ret == mbox.Yes:
    self.centralLogView.launchCmdIntoPopen(cmd)
    return True


  ##################################################################
  def getStrLaunchBat(self, etudeDir, args="4bit etc."):
    # data Corteo in etudeDir/../data
    # logger.debug("current etude directory %s" % etudeDir)

    code_exe = os.path.join(*"${SOLVERLABGUI_ROOT_DIR} solverlabCode solverlab_mingw64.exe".split())
    code_exe = os.path.realpath(os.path.expandvars(code_exe))
    data_dir = os.path.join(etudeDir, "..", "data")
    data_dir = os.path.realpath(os.path.expandvars(data_dir))

    logger.info("solverlab code file %s" % code_exe)
    logger.info("solverlab data dir  %s" % data_dir)

    if not os.path.isdir(data_dir):
      msg = "rem inexisting corteo data directory: %s\nrem fix it." % data_dir
      logger.critical(msg)
      return msg

    mater_file = "Materials.in"
    struc_file = "Structure.in"
    compo_file = "Composition.in"

    material_file = os.path.join(etudeDir, mater_file)
    structure_file = os.path.join(etudeDir, struc_file)
    composition_file = os.path.join(etudeDir, compo_file)
    cmd = """\
@echo off
SET studyDir={0}
SET solverlabExe={1}
CD %studyDir%

REM tree /A
REM dir

REM to see solverlab progress (if log as file)
REM type {0}/solverlab.log

REM get solverlab log as file
REM %solverlabExe% -p 9 -data ../data -c ./Configuration.in > ./solverlab.log

REM get solverlab log as stdout
%solverlabExe% -p 9 -data ../data -c ./Configuration.in
SET solverlabExitCode=%errorlevel%

ECHO END of solverlab, exit code is %solverlabExitCode%
EXIT %solverlabExitCode%
""".format(join(etudeDir), join(code_exe))
    return cmd

  ##################################################################
  def getStrLaunchSh(self, etudeDir, args="4bit etc."):
    # data Corteo in etudeDir/../data
    # logger.debug("current etude directory %s" % etudeDir)

    PACKAGESPY_ROOT_DIR = os.getenv("PACKAGESPY_ROOT_DIR")
    SOLVERLABGUI_ROOT_DIR = os.getenv("SOLVERLABGUI_ROOT_DIR")
    nbproc = int(self._model.Analysis.caseSolverlab.NumberOfProcessors)
    if nbproc == 0:
      nbproc = 1
    """
    if PACKAGESPY_ROOT_DIR is not None: # as salome_matix installation
      code_exe = os.path.join(*"${PACKAGESPY_ROOT_DIR} bin solverlab_cea.exe".split())
    else: # as default
      code_exe = os.path.join(*"${SOLVERLABGUI_ROOT_DIR} solverlabCode solverlab_linux64.exe".split())

    code_exe = os.path.realpath(os.path.expandvars(code_exe))

    data_dir = os.path.join(etudeDir, "..", "data")
    data_dir = os.path.realpath(os.path.expandvars(data_dir))

    logger.info("solverlab code file %s" % code_exe)
    logger.info("solverlab data dir  %s" % data_dir)

    if not os.path.isdir(data_dir):
      msg = "# inexisting corteo data directory: %s\n# fix it." % data_dir
      logger.critical(msg)
      return msg
    """
    cmd = f"""
#!/usr/bin/env bash

# This script is generated automatically.
# PACKAGESPY_ROOT_DIR environ variable should be set before use
# $SOLVERLAB_INSTALL/env_solverlab.sh

echo "########## START"

export PACKAGESPY_ROOT_DIR={PACKAGESPY_ROOT_DIR}
export SOLVERLABGUI_ROOT_DIR={SOLVERLABGUI_ROOT_DIR}

cd {etudeDir}

if [ -f lock ]; then
    echo "lock file exist and Solverlab supposed working in directory."
    echo "delete file lock if you are sure."
    cat lock
    echo "ERROR on locked Solverlab directory"
    echo "########## END"
    exit 2
fi

echo "hostname="`hostname` > lock
echo "launchTime="`date` >> lock
echo "launchDirectory="`pwd` >> lock
rm -f lock_old

python {SOLVERLABGUI_ROOT_DIR}/solverlabpy/runSolverlab.py solverlabGui.xml

if [ $? -ne 0  ]; then
    echo "run test failed !"
fi

echo "pidSolverlab="$! >> lock
cat lock
wait

echo "endTime="`date` >> lock
mv -f lock lock_old
echo "Solverlab finished"

echo "########## END"

"""
    return cmd

  def SaveSolverlabModelXmlAction(self, withMessage=True):
    if self._model == None:
      QTW.QMessageBox.warning(self._desktop, "warning", "No Solverlab data to save")
      return False
    etudeDir = self._model.getEtudeWorkdirExpanded()
    logger.debug("SaveSolverlabModelXmlAction %s" % etudeDir)
    if not os.path.isdir(etudeDir):  # create without question if inexisting
      self.CreateEtudeWorkdir()
    nameFile = os.path.join(etudeDir, self._solverlabXmlDefaultName)
    self._model.toFileXml(nameFile)
    if withMessage:
      QTW.QMessageBox.information(self._desktop, "info", "saved file:\n'%s'" % nameFile)
    return True

  def LoadSolverlabDefaultModelAction(self):
    logger.debug("LoadSolverlabDefaultModel")
    if True: #try:
      aData = ModelSvl()
      aData.setDefaultValues()
      self.clearModel()
      self.setModel(aData)
      self.RefreshSolverlabModelSignal.emit(None)
      #self._model.userExpand()
      self.ExpandAllAction()
      return True
    else: #except:
      import traceback
      trace = traceback.format_exc()
      # traceback.print_exc() #better explicit verbose problem
      QTW.QMessageBox.warning(self._desktop, "warning",
                                "Load Solverlab Default data: houston! we have a problem \n\n%s" % trace)
      return False

  def LoadSolverlabModelXmlAction(self):
    from widgetpy.salomeQFileDialog import SalomeQFileDialog
    aDialog = SalomeQFileDialog(parent=self._desktop)
    aDir = USET.getExpandedVar("_SOLVERLABGUI_WORKDIR")
    nameFile = aDialog.browseFileDialog('Load Solverlab file xml', aDir, "(*.xml *.XML)", [])
    logger.debug("LoadSolverlabModelXml %s" % nameFile)
    if nameFile == "": return True  # cancel
    realPath = os.path.realpath(nameFile)
    if not os.path.isfile(realPath):
      QTW.QMessageBox.warning(self._desktop, "warning",
                                "Load Solverlab xml data: not a file \n'%s'" % realPath)
      return False
    try:
      aData = UXYZ.fromFileXml(realPath)
      logger.debug("LoadSolverlabModelXml: data loaded from file")
      if aData.__class__.__name__ != "ModelSvl":
        QTW.QMessageBox.warning(self._desktop, \
          "warning", "Load Solverlab data: '%s' is not a ModelSvl instance from \n'%s'" % (aData.__class__.__name__, realPath))
        return False
      self.clearModel()
      self.setModel(aData)
      # logger.debug("LoadSolverlabModelXml: set model data done")
      self.RefreshSolverlabModelSignal.emit(None)
      self._model.userExpand()
      self.ExpandAllAction()
      return True
    except Exception as e:
      # traceback.print_exc()
      QTW.QMessageBox.warning(self._desktop, "warning","Load Solverlab xml data: problem loading\n%s" % e)
      return False

  def loadXmlFile(self, aFileOrDir):
    realPath = os.path.expandvars(aFileOrDir)
    if os.path.isdir(realPath):
      realPath = os.path.join(realPath, "solverlabGui.xml")  # default name
    if not os.path.isfile(realPath):
      logger.error("inexisting file: %s" % realPath)
      return

    etudeDirToLoad, _ = os.path.split(realPath)
    USET.setEnvVar("SOLVERLABGUIDIR", etudeDirToLoad)
    aData = UXYZ.fromFileXml(realPath)

    self.clearModel()
    self.setModel(aData)
    self.checkFromCpack()
    self.RefreshSolverlabModelSignal.emit(None)
    self.ExpandAllAction()

  def isPresentInFile(self, aFile, strOrListOfStr):
    """return True if str or one str of [str1, str2,...] is in file"""
    try:
      with open(aFile, "r") as f:
        contents = f.read()
      if type(strOrListOfStr) == str:
        if strOrListOfStr in contents: return True
      elif type(strOrListOfStr) == list:
        for aStr in strOrListOfStr:
          if aStr in contents: return True
    except:
      pass
    return False

  def _UserExpand(self, aList):
    logger.debug("User expand\n%s" % PP.pformat(aList))
    for treeView in self.treeViews:
      treeView.userExpandModel(aList)
    return True

  def _UserMode(self, aStr):
    import solverlabpy.configSvl as CFGSVL
    logger.debug("User change mode '%s'" % aStr)
    CFGSVL.setCurrentMode(aStr)
    self.RefreshSolverlabModelSignal.emit(None)
    return True


def launchFromSalomePyConsole():
  from solverlabpy.controllerSvl import ControllerSvl
  desktop = QTW.QMainWindow()
  desktop.resize(1000, 700)
  ctrl = ControllerSvl(desktop=desktop)
  desktop.show()
  return desktop


class OutWrapper(object):
  def __init__(self, edit):
    self.out = sys.stdout
    self.textEdit = edit

  def write(self, message):
    self.out.write(message)
    self.textEdit.insertPlainText(message)

  def flush(self):
    self.out.flush()


if __name__ == '__main__':
  from salomepy.onceQApplication import OnceQApplication

  app = OnceQApplication()
  aWindow = launchFromSalomePyConsole()
  app.exec_()
