# -*- coding: utf-8 -*-
#  Copyright (C) 2007-2010  CEA/DEN, EDF R&D, OPEN CASCADE
#
#  Copyright (C) 2003-2007  OPEN CASCADE, EADS/CCR, LIP6, CEA/DEN,
#  CEDRAT, EDF R&D, LEG, PRINCIPIA R&D, BUREAU VERITAS
#
#  This library is free software; you can redistribute it and/or
#  modify it under the terms of the GNU Lesser General Public
#  License as published by the Free Software Foundation; either
#  version 2.1 of the License.
#
#  This library is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#  Lesser General Public License for more details.
#
#  You should have received a copy of the GNU Lesser General Public
#  License along with this library; if not, write to the Free Software
#  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA
#
#  See http://www.salome-platform.org/ or email : webmaster.salome@opencascade.com
#
# Author : A. Bruneton
#

from PyQt4.QtCore import Qt 
from PyQt4.QtGui import QMainWindow,QMenu, QDockWidget
from MainCFWidget import MainCFWidget

class CFDesktop(QMainWindow):
    """ 
    """
    VIEW_TYPE = "COREFLOWS"
        
    def __init__(self, sgPyQt):
        QMainWindow.__init__(self)
        self._sgPyQt = sgPyQt
        self._mainView = None
        self._viewID = -1
        
    def initialize(self):
        """ Initialize is called later than __init__() so that the Desktop and the SgPyQt
        objects can be properly initialized.
        """
        self._currID = 1235
        
        self._sgDesktop = self._sgPyQt.getDesktop()
        self.createIDs()
        self.createActions()
        
        self.createToolbars()
        self.createMenus()

    def generateID(self):
        self._currID += 1
        return self._currID
    
    def createIDs(self):
        pass
        # Actions
#        self.itemDelActionID = self.generateID()
#        self.cpsActionID = self.generateID()
#        self.addPSActionID = self.generateID()
#        
#        # Menus
#        self.etudeMenuID = self.generateID()

    def createActions(self):
      pass
#        ca = self._sgPyQt.createAction
#        self.itemDelAction = ca(self.itemDelActionID, "Delete selected", "Delete selected", "", "")
#        self.cpsAction = ca(self.cpsActionID, "Clear plot set", "Clear plot set", "", "")
#        self.addPSAction = ca(self.addPSActionID, "Add plot set", "Add plot set", "", "")

    def createToolbars(self):
        pass
#         self.Toolbar = self._sgPyQt.createTool(self.tr("Toolbar"))
#         self._sgPyQt.createTool(self.fileNewAction, self.Toolbar)
#         self._sgPyQt.createTool(self.filePrintAction, self.Toolbar)
#         sep = self._sgPyQt.createSeparator()
#         self._sgPyQt.createTool(sep, self.Toolbar)
#         self._sgPyQt.createTool(self.editUndoAction, self.Toolbar)
#         self._sgPyQt.createTool(self.editRedoAction, self.Toolbar)

    def createMenus(self):
      pass
#      curveMenu = self._sgPyQt.createMenu( "Plot management", -1, self.etudeMenuID, self._sgPyQt.defaultMenuGroup() )
#      self._sgPyQt.createMenu(self.itemDelAction, curveMenu)

    def createView(self):
      if self._mainView is None:
        self._mainView = MainCFWidget()
      vid = self._sgPyQt.createView(self.VIEW_TYPE, self._mainView)
      return vid
          
    def showCentralWidget(self):
      if self._viewID == -1:
        self._viewID = self.createView()
      else:
        self._sgPyQt.activateView(self._viewID)
        
    def hideCentralWidget(self):
      if self._viewID != -1:
        self._sgPyQt.setViewVisible(self._viewID, False)
