#!/usr/bin/env python
# -*- coding: utf-8 -*-

# %% LICENSE_SALOME_CEA_BEGIN
# see PACKAGESPY/LICENCE file
# %% LICENSE_END

import os

from PyQt5.QtGui import *
from PyQt5.QtCore import *
from PyQt5.QtWidgets import *

def verbose():
  return True

"""class sgPyQtEmulate( QObject ):
  def __init__(self):
    QObject.__init__(self)
    pass"""
    
class DesktopEmulateSalome( QMainWindow ) :        
  def __init__( self, title="" ):     
    QMainWindow.__init__(self)
    self.setWindowTitle( title+"Salome" )
    self.setMinimumWidth( 800 )
    self.setMinimumHeight( 600 )
    self.myActions = {}
    self.myMenus = {}
    self.menubar = self.menuBar()
    self.statusbar = self.statusBar()
    self.dict_command = {}

  def getDesktop(self):
    return self
  
  def createAction( self, idSalome, title, tooltip, statustip, image=None ):
    if verbose(): print('create action', idSalome, title, image)
    if image == None:
      anAction=QAction(title, self)
    else:
      icon = QIcon( os.path.join( "..", "..", "resources", image ) )
      anAction=QAction(icon, title, self)
    anAction.setObjectName( title )
    anAction.setToolTip( tooltip )
    anAction.setStatusTip( statustip )
    self.addAction( anAction )
    self.myActions[idSalome] = anAction
    if verbose(): print("***DesktopEmulateSalome.dict_command",list(self.dict_command.keys()))
    anAction.triggered.connect( self.dict_command[idSalome] )
    return anAction

  def createMenu( self, p1, p2, p3=None, p4=None ):
    """
    (title, noMenu, idSalome, menugroup) or
    (action, menu)
    """
    if p3==None:
      action, menu = p1, p2
      if verbose(): print('create action', action.objectName(), 'in menu', menu.objectName())
      menu.addAction(action)  
      return None
    else:  
      title, noMenu, idSalome, menugroup = p1, p2, p3, p4
      if verbose(): print('create menu vierge', title, noMenu, idSalome, menugroup)
      self.myMenus[noMenu] = self.menubar.addMenu(title)
      self.myMenus[noMenu].setObjectName( title )
      return self.myMenus[noMenu]
          
  def createTool( self, p1, p2=None):
    """
    (title) or
    (action, toolbar)
    """
    if p2==None:
      title = p1
      if verbose(): print('create toolbar', title)
      toolbar = self.addToolBar(title)  
      toolbar.setObjectName( title )
      return toolbar
    else:  
      action, toolbar = p1, p2
      if verbose(): print('add action', action.objectName(), 'in toolbar', toolbar.objectName())
      toolbar.addAction(action)
      return None
          
  def defaultMenuGroup(self):
    return self.menubar
  
  def stringSetting( self, nameModule, what, defaultWhat):
    return defaultWhat
