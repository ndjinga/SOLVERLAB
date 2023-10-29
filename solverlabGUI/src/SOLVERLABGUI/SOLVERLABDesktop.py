#!/usr/bin/env python
# -*- coding: utf-8 -*-

# %% LICENSE_SALOME_CEA_BEGIN
# see SOLVERLABGUI/LICENSE file
# %% LICENSE_END


import os
from PyQt5.QtGui import *
from PyQt5.QtCore import *
from PyQt5.QtWidgets import *

import SOLVERLAB_utils as MUT # as Module_UTils
from solverlabpy.controllerSvl import ControllerSvl as Controller

_NAME = "SOLVERLAB"
_NAMEGUI = "SOLVERLABGUI"
_Name = _NAME.title()

def verbose():
  return True

def log(mess):
  if verbose: MUT.log_print(_NAME+ "Desktop" + mess)

class SOLVERLABDesktop( QObject ):

   _ORIGIN_ID = MUT.moduleID() + 100
   _SIMULATION_MENU_ID = _ORIGIN_ID + 0
   #some here and do not pass _END_ID
   #_NEW_SIMULATION_ID = _ORIGIN_ID + 1
   _END_ID =  _ORIGIN_ID + 10

   def __init__( self, sgPyQt, sg  ):

       QObject.__init__( self )
       self._sgPyQt = sgPyQt
       self._sg = sg #swig salome
       self._sgDesktop = self._sgPyQt.getDesktop()
       self._widgetDialogBox = None
       self._controller = Controller( self )
       self.createTreeView()
       self.addTreeDockDone = False
       self._globalSimplePlotView = None
       self._globalWindowForLog = None
       self._globalSimplePlotViewID = None
       self._globalWindowForLogID = None
       #Getting the help file
       self._helpFile = os.path.join( os.getenv(_NAMEGUI + "_ROOT_DIR"),"doc","index.html" )
       self._actions = []
       self.config = None  # user configuration will come from file ${SOLVERLAB_WORKDIR}/*cfg

   def tryHide(self, anObject):
       """avoid RuntimeError: wrapped C/C++ object of type QDockWidget has been deleted"""
       if anObject == None:
         return False
       try:
         anObject.hide()
         return True
       except:
         return False

   def tryShow(self, anObject):
       """avoid RuntimeError: wrapped C/C++ object of type QDockWidget has been deleted"""
       if anObject == None:
         return False
       try:
         anObject.show()
         return True
       except:
         return False

   def set_dict_command( self, dict_command ):
       ii = self._END_ID  #marging if user want to set some actions directly
       if self._controller == None: return
       for action in self._controller.actions:
         # log("Add controller action '%s' in dict_command of salome desktop" % action.name)
         theSlot = action.slot
         dict_command[ ii ] = theSlot
         ii += 1

   def findDockByName(self, name):
       children = self._sgDesktop.children()
       # log("findDockByName %s" % name)
       for c in children:
         if isinstance(c, QDockWidget):
           # log(" windowTitle() %s %s" % (c.windowTitle(), id(c)))
           if c.windowTitle() == name: return c
       log(" findDockByName: dock '%s' not found" % name)
       return None

   def tabifyDockWidgets(self, toTabify):
       """
       toTabify is a list of titles of QDockWidgets or the QDockWidgets themselves
       lenght have to be > 2 dockwidgets toTabify
       """
       children = self._sgDesktop.children()
       if len(toTabify) < 2:
         log(" tabifyDockWidgets: 2 docks or more")
         return
       aFirstDock = toTabify[0]
       if isinstance(aFirstDock, str): aFirstDock = self.findDockByName(aFirstDock)
       if aFirstDock == None:
         log(" tabifyDockWidgets: first Dock to tabify not found:",toTabify[0])
         #on pourrait prendre le premier trouve? non.
         return

       for aDock in toTabify[1:]:
         if isinstance(aDock, str):
           theDock = self.findDockByName(aDock)
         else:
           if isinstance(aDock, QDockWidget):
             theDock = aDock
           else:
             theDock = None
         if theDock == None : continue
         self._sgDesktop.tabifyDockWidget(aFirstDock, theDock)

   def getWidgetDialogBox(self):
       #only one, which change title and inside
       if self._widgetDialogBox == None:
         self._widgetDialogBox = QDockWidget( self._sgDesktop )
       return self._widgetDialogBox


   def createTreeView( self, treeView=None ):
       log(" createTreeView")
       self._globalTree = self._controller.treeViews[0] #TreeWidget( self )
       self._dockGlobalTree = self._controller.docks[0]
       """
       self._globalTree.header().close()
       self._dockGlobalTree = QDockWidget( _Name + "Objects", self._sgDesktop )
       self._dockGlobalTree.setWidget( self._globalTree )
       """

   def activateDesktop( self ):
       log(" activateDesktop done: %s" % self.addTreeDockDone)
       import solverlabpy.configSvl as CFGSVL
       self.config = CFGSVL.getMainConfig()
       log('config\n%s' % self.config )

       if self.addTreeDockDone is False:
         self.createActions()
         self.createMenus()
         self.createToolBars()
         self.createPopups()
         res = self._sgDesktop.addDockWidget( Qt.LeftDockWidgetArea, self._dockGlobalTree )
         log(" activateDesktop addDockWidget %s" % res)
         self._sgDesktop.setTabPosition(Qt.LeftDockWidgetArea, QTabWidget.North)
         self.addTreeDockDone=True
         #self.createPlotView()
         #self.createWindowForLog()

       self.getGlobalSimplePlotView()
       self.getGlobalWindowForLog()
       if self._globalSimplePlotViewID is not None:
         self._sgPyQt.setViewVisible(self._globalSimplePlotViewID, True)
       if self._globalWindowForLogID is not None:
         self._sgPyQt.setViewVisible(self._globalWindowForLogID, True)

       self._dockGlobalTree.show()
       self.tabifyDockWidgets(["Object Browser", _Name + "Objects"])
       toClose = self.findDockByName("Object Browser")
       self.tryHide(toClose)
       #toShow = self.findDockByName("Object Browser")
       #self.tryShow(toShow)
       toClose = self.findDockByName("NoteBook")
       self.tryHide(toClose)
       toShow = self.findDockByName(_Name + "Objects")
       self.tryShow(toShow)
       if self._globalWindowForLogID != None:
         self._sgPyQt.setViewVisible(self._globalWindowForLogID, True)

   def deactivateDesktop( self ):
       log(" deactivateDesktop")
       toShow = self.findDockByName("Object Browser")
       self.tryShow(toShow)
       self.tryHide(self._dockGlobalTree) #may be user delete
       self._sgPyQt.setViewVisible(self._globalSimplePlotViewID, False)
       self._sgPyQt.setViewVisible(self._globalWindowForLogID, False)
       #self.tryHide(self._globalSimplePlotView)
       #self.tryHide(self._globalWindowForLog)

   def createGraphicsView( self ):
       scene = GraphicsScene( self._controller )
       self._globalGraphicsView = GraphicsView( scene )
       self._globalGraphicsViewID = self._sgPyQt.createView( "ViewCurve", self._globalGraphicsView )

   def getGlobalSimplePlotView(self):
       try:
         self._globalSimplePlotView.on_draw()
       except:
         #sometimes gui user delete viewer, recreate it!
         self.createPlotView()
       return self._globalSimplePlotView

   def getGlobalWindowForLog(self):
       if self._globalWindowForLog is None:
         self.createWindowForLog()

       try:
         self._globalWindowForLog.show()
         ok = True
         return self._globalWindowForLog
       except Exception as e:  # RuntimeError: wrapped C/C++ object of type QMainWindowForLog has been deleted
         # sometimes gui user delete viewer, recreate it!
         ok = False
         log("ERROR: getGlobalWindowForLog problem\n%s" % e)

       if not ok: # sometimes gui user delete viewer, recreate it!
         self._globalWindowForLog = None
         self._controller.centralLogView = None
         self.createWindowForLog()
         self._globalWindowForLog.show()
         return self._globalWindowForLog

   def createPlotView( self ):
       from matplotlibpy.matplotlibWindowToolbar import MatplotlibWindowToolbar
       self._globalSimplePlotView = MatplotlibWindowToolbar()
       self._globalSimplePlotViewID = self._sgPyQt.createView( _Name + "Plot", self._globalSimplePlotView )
       log(" createPlotView controller: %i view: %i" % ( id(self._controller), self._globalSimplePlotViewID ))

   def createWindowForLog( self ):
       if self._controller.centralLogView is None: # sometimes gui user delete viewer, recreate it!
         self._controller.addCentral()
       centralLogView = self._controller.centralLogView
       for dock in centralLogView.docks:
         self.tryHide(dock) # only tabs, not treeview
       self._globalWindowForLog = centralLogView
       self._globalWindowForLogID = self._sgPyQt.createView( _Name + "Log", self._globalWindowForLog )
       log(" createWindowForLog controller: %i view: %i" % ( id(self._controller), self._globalWindowForLogID ))

   def activatePlotView( self ):
       if self._globalSimplePlotViewID != None:
         self._sgPyQt.activateView( self._globalSimplePlotViewID)

   def activateLogView( self ):
       if self._globalWindowForLogID != None:
         self._sgPyQt.activateView( self._globalWindowForLogID )

   def createActions( self ):
       ii = self._END_ID  #marging if user want to set some actions directly
       if self._controller == None: return
       for action in self._controller.actions:
         # log(" Add controller action '%s' in salome desktop" % action.name)
         self._actions.append(
           self._sgPyQt.createAction( ii, action.name, action.text(), action.toolTip(), action.iconFile ) )
         ii += 1

   def createMenus( self ):
       simulationMenu = self._sgPyQt.createMenu( _Name, -1, self._SIMULATION_MENU_ID, self._sgPyQt.defaultMenuGroup() )
       ii = self._END_ID  #marging if user want to set some actions directly
       if self._controller == None: return
       for action in self._controller.actions:
         # log(" Add controller action '%s' in salome desktop menu" % action.name)
         self._sgPyQt.createMenu( action, simulationMenu )
         ii += 1

   def createToolBars( self ):
       simulationTB = self._sgPyQt.createTool(_Name)
       ii = self._END_ID  #marging if user want to set some actions directly
       if self._controller == None: return
       for action in self._controller.actions:
         # log(" Add controller action '%s' in salome desktop toolbar" % action.name)
         self._sgPyQt.createTool( action, simulationTB )
         ii += 1

   def createPopups( self ):
       pass

   def getController( self ):
       return self._controller

   def setController( self, controller ):
       self._controller = controller

   def getGlobalTree( self ):
       return self._globalTree

   def getDockGlobalTree( self ):
       return self._dockGlobalTree
