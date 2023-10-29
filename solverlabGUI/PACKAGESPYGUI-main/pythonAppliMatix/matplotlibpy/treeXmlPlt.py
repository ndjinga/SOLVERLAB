#!/usr/bin/env python
# -*- coding: utf-8 -*-

# %% LICENSE_SALOME_CEA_BEGIN
# see PACKAGESPY/LICENCE file
# %% LICENSE_END


from PyQt5 import QtCore, QtGui, QtWidgets
from xyzpy.guiXyz.treeXmlXyz import TreeXmlXyz
from salomepy.browserPythonClass import BrowserPythonClass

verbose = False
verboseEvent = False

"""
cosmetic stuff for actions on plots Figures matplotlib (plt)
"""

class TreeXmlPlt(TreeXmlXyz):
  
  class COLS:
    labels = ['Name', 'Value', 'Attributes']
    Tag = 0
    Text = 1
    Attributes = 2

  def __init__(self, parent=None):
    super(TreeXmlPlt, self).__init__(parent)
    
    self.setHeaderLabels(self.COLS.labels)
    font = QtGui.QFont(self.font()) #copy: self.font() is const
    font.setFamily("Monospace")
    font.setPointSize(9)
    self.setFont(font)
    self.setAlternatingRowColors(True)
    pal=self.palette()
    pal.setColor(pal.Base, QtGui.QColor(210,230,220))
    pal.setColor(pal.Text, QtGui.QColor(0,0,0))
    self.setPalette(pal)
    
  def _createActions(self):
    super(TreeXmlPlt, self)._createActions()
    self.pltActions = [
      self._createAction('Display plot', None, 'Display plot', self.displayPlt, None),
      self._createAction('Display standalone plot', None, 'Display standalone plot', self.displayStandalonePlt, None),
      #self._createAction('Delete plot', None, 'Delete plot', self.deletePlt, None),
      self._createAction('Delete plot', None, 'Delete plot', self.deletePlt, None),
      ]

  def _createContextMenus(self):
    menus = super(TreeXmlPlt, self)._createContextMenus()
    menuPlt = QtWidgets.QMenu("Plot", self)
    for action in self.pltActions: 
      menuPlt.addAction(action)
      menus[0].addAction(action)
    menus[1] = menuPlt
    return menus
    
  def _getNameOrIdentFigureFromItem(self,  item):
    text0 = str(item.text(0)).strip()
    if text0 == "name": 
      return str(item.text(1))
    if text0 == "ident": 
      return str(item.text(1))
    else:
      return text0

  def displayStandalonePlt(self):
    items = self.selectedItems()
    if len(items) == 0:
      QtWidgets.QMessageBox.warning(self, "warning", "select one plot")
      return
    controller = self.getController()
    #self._tmp = []
    for item in items:
      #tmp = BrowserPythonClass(item, levelMax=1)
      #tmp.display()
      #self._tmp.append(tmp)
      #print 'display item selected', item.text(0), type(item) #item.text(1), item.text(2)
      ok,  why = controller.displayStandalonePlt(self._getNameOrIdentFigureFromItem(item))
      if not ok:
        QtWidgets.QMessageBox.warning(self, "warning", why)
      return

  def displayPlt(self):
    items = self.selectedItems()
    if len(items) == 0:
      QtWidgets.QMessageBox.warning(self, "warning", "select one plot")
      return
    if len(items) > 1:
      QtWidgets.QMessageBox.warning(self, "warning", "select only one plot")
      return
    controller = self.getController()
    for item in items:
      ok,  why = controller.displayPlt(self._getNameOrIdentFigureFromItem(item))
      if not ok:
        QtWidgets.QMessageBox.warning(self, "warning", why)
      return

  def deletePlt(self):
    items = self.selectedItems()
    if len(items) == 0:
      QtWidgets.QMessageBox.warning(self, "warning", "select one plot")
      return
    controller = self.getController()
    for item in items:
      #print 'delete item selected', item.text(0), type(item) #item.text(1), item.text(2)
      ok,  why = controller.deletePlt(self._getNameOrIdentFigureFromItem(item))
      if not ok:
        QtWidgets.QMessageBox.warning(self, "warning", why)
      return


