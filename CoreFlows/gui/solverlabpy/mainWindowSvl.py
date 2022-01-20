#!/usr/bin/env python
# -*- coding: utf-8 -*-

# %% LICENSE_SALOME_CEA_BEGIN
# Copyright (C) 2008-2018  CEA/DEN
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


import sys
import pprint as PP
import json

from PyQt5 import QtGui, QtCore, QtWidgets as QTW
import ribbonpy.ribbonIcons as IUSR
from ribbonpy.qTabMultipleText import QTabMultipleText
from ribbonpy.ribbonWidget import RibbonWidget
from ribbonpy.ribbonQMainWindow import QMainWindowForRibbon

verbose = True
verboseEvent = False

from ribbonpy.ribbonTrace import getLoggerRibbon
RT = getLoggerRibbon()


class QMainWindowSvl(QMainWindowForRibbon):
  """
  QMainWindow with ribbonWidget for solverlab
  """

  index = [0]  # for unique instance naming, if some

  def __init__(self, *args, **kwargs):
    super(QMainWindowSvl, self).__init__(*args)
    self.setObjectName("QMainWindowSvl%s" % str(self.index))
    self.setWindowTitle("SOLVERLAB GUI") # default

  def setFromJson(self, valuesJson=None):
    """
    virtual, only for example
    values not used, yet
    """
    RT.info("QMainWindowForSolverlab.setFromJson: create widget")
    self.__addCentral()
    self.__createActions()
    self.__addToolBars()
    # self.__addRibbonDock()
    # self.__addRibbonToolBar()
    # self.__addDocks()

    self.statusBar().showMessage('Ready')
    # self.statusBar().addPermanentWidget(QTW.QLabel("QLabel (or else...) for example"), stretch=0)
    # self.resize(900, 500)

  """
  def close(self):
    #warning: is not usable if salome, only catch hide event, not close
    print "QMainWindowForSolverlab %s close",self.objectName()
    return super(QMainWindowForSolverlab, self).close()

  def closeEvent(self, event):
    print "QMainWindowForSolverlab %s closeEvent",self.objectName()
    #event.ignore()
    return super(QMainWindowForSolverlab, self).closeEvent(event)

  def event(self, event):
    print "QMainWindowForSolverlab %s event",self.objectName(),strEvent(event)
    return super(QMainWindowForSolverlab, self).event(event)

  def __del__(self):
    print "QMainWindowForSolverlab %s __del__",self.objectName()
    return super(QMainWindowForSolverlab, self).__del__()
  """

  def __addCentral(self, valuesJson=None):
    """
    virtual, only for example
    values not used, yet
    """
    RT.info("QMainWindowForSolverlab.__addCentral: create centralWidget QTabMultipleText")
    self.central = QTabMultipleText()
    self.centralWidget()
    self.setCentralWidget(self.central)
    # self.centralWidget().resize(self.centralWidget().size())
    self.centralWidget().show()

  def __addDocks(self, valuesJson=None):
    """
    virtual, only for example
    values not used, yet
    """
    RT.info("QMainWindowForSolverlab.__addDocks: For example only")
    self.docks = []
    dock = QTW.QDockWidget("TreeForExample", self)
    self.treeView = QTW.QTreeWidget()
    dock.setWidget(self.treeView)
    pos = QtCore.Qt.LeftDockWidgetArea
    dock.setAllowedAreas(pos)
    self.docks.append(dock)
    for dock in self.docks:
      self.addDockWidget(pos, dock)

  def __addRibbonDock(self, valuesJson=None):
    """
    virtual, only for example
    values not used, yet
    """
    RT.info("QMainWindowForSolverlab.__addRibbonDock: create Ribbon in dock")
    import ribbonpy.ribbonClassFactory as RCF
    self.ribbon = RibbonWidget(setFromJson=RCF.getExampleJsonRibbon())
    dock = QTW.QDockWidget(self.ribbon.objectName() + "_dock", self)
    p = dock
    # p.setWindowFlags(QtCore.Qt.CustomizeWindowHint)
    p.setWindowFlags(QtCore.Qt.FramelessWindowHint)
    # p.setFeatures(p.DockWidgetVerticalTitleBar)
    # p.setWindowFlags(QtCore.Qt.Tool)
    dock.setWidget(self.ribbon)
    pos = QtCore.Qt.TopDockWidgetArea
    dock.setAllowedAreas(pos)
    self.docks.append(dock)
    for dock in self.docks:
      self.addDockWidget(pos, dock)
      dock.close()

  def __addRibbonToolBar(self, valuesJson=None):
    """
    virtual, only for example
    values not used, yet
    """
    RT.info("QMainWindowForSolverlab.__addRibbonToolBar: create Ribbon in ToolBar")
    import ribbonpy.ribbonClassFactory as RCF
    self.ribbon2 = RibbonWidget(setFromJson=RCF.getExampleJsonRibbon())
    tb = self.addToolBar("ribbon2_toolbar")
    tb.addWidget(self.ribbon2)
    tb.setAllowedAreas(QtCore.Qt.TopToolBarArea | QtCore.Qt.BottomToolBarArea)
    self.toolBars.append(tb)

  def __addToolBars(self, valuesJson=None):
    """
    virtual, only for example
    values not used, yet
    """
    RT.info("QMainWindowForSolverlab.__addToolBars: not yet implemented")
    return
    self.toolBars = []
    tb = self.addToolBar("toolbar1")
    for action in self.actions:
      tb.addAction(action)
    self.toolBars.append(tb)

  def __createActions(self, valuesJson=None):
    """
    virtual, only for example
    values not used, yet
    """
    """create general actions"""
    RT.info("QMainWindowForSolverlab.__createActions: not yet implemented")
    return
    self.actions = []
    self.actions.append(self.__createAction("Display File", "D", "Display file", self.central.displayFile, "open"))

  def __createAction(self, Name, Shortcut, ToolTip, Call, Icon=None):
    """create one actions"""
    # http://pyqt.sourceforge.net/Docs/PyQt5/signals_slots.html
    action = QTW.QAction(Name, self)
    if Shortcut != None: action.setShortcut(self.prefixShortcut + Shortcut)
    action.setToolTip(ToolTip)
    if Icon != None:
      action.setIcon(IUSR.getIconFromName(Icon))
    action.triggered.connect(Call)
    return action


