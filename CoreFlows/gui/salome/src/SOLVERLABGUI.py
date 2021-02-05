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
# Author : A. Bruneton and B. Secher
#
import traceback
import string
import os
import sys
from qtsalome import *

import salome
from CFDesktop import CFDesktop

# Get SALOME PyQt interface
import SalomePyQt
import libSALOME_Swig

########################################################
# Global variables
########################################################

sgPyQt = SalomePyQt.SalomePyQt()
sg = libSALOME_Swig.SALOMEGUI_Swig()
sgDesktop = sgPyQt.getDesktop()
widgetDialogBox = None

moduleDesktop  = None

########################################################
# Internal methods
########################################################

def getDesktop():
    """This method returns the current SOLVERLABT desktop"""

    global moduleDesktop
    return moduleDesktop

def setDesktop( ):
    """This method sets and returns SOLVERLABT desktop"""

    global moduleDesktop

    if moduleDesktop is None:
        moduleDesktop = CFDesktop( sgPyQt )
        pass
    return moduleDesktop

def incObjToMap( m, id ):
    """This method incrementes the object counter in the map"""

    if id not in m: m[id] = 0
    m[id] += 1
    pass

def getSelection():
    """This method analyses selection"""

    selcount = sg.SelectedCount()
    seltypes = {}
    for i in range( selcount ):
        incObjToMap( seltypes, getObjectID( sg.getSelected( i ) ) )
        pass
    return selcount, seltypes

################################################
# Callback functions
################################################

def initialize():
    """This method is called when module is initialized. It performs initialization actions"""
    setDesktop()
    pass

def windows():
    """This method is called when module is initialized. It returns a map of popup windows to be used by the module"""

    wm = {}
    wm[SalomePyQt.WT_ObjectBrowser] = Qt.LeftDockWidgetArea
    wm[SalomePyQt.WT_PyConsole]     = Qt.BottomDockWidgetArea
    return wm

def views():
    """This method is called when module is initialized. It returns a list of 2D/3D views to be used by the module"""
    return []

def createPreferences():
    """This method is called when module is initialized. It exports module preferences"""
    pass

def activate():
    """This method is called when module is initialized. It returns True if activating is successfull, False otherwise"""

    global moduleDesktop

    fv = moduleDesktop.showCentralWidget()
    return True

def viewTryClose( wid ):
    sgPyQt.setViewClosable(wid, True)
    pass

def deactivate():
    """This method is called when module is deactivated"""

    global moduleDesktop, widgetDialogBox
    moduleDesktop.hideCentralWidget()
    pass

def activeStudyChanged():
    """This method is called when active study is changed"""

    setDesktop()
    pass

def createPopupMenu( popup, context ):
    """This method is called when popup menu is invocked"""
    pass

def OnGUIEvent( commandID ):
    """This method is called when a GUI action is activated"""

    if commandID in dict_command:
       dict_command[commandID]()
       pass
    pass

def preferenceChanged( section, setting ):
    """This method is called when module's preferences are changed"""
    pass

def activeViewChanged( viewID ):
    """This method is called when active view is changed"""
    pass

def viewCloned( viewID ):
    """This method is called when active view is cloned"""
    pass

def viewClosed( viewID ):
    """This method is called when active view viewClosed"""
    pass

def engineIOR():
    """This method is called when study is opened. It returns engine IOR"""
    return getEngineIOR()

########################################################
