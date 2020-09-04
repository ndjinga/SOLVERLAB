# -*- coding: latin-1 -*-
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
import sys, os
from CFDesktop import CFDesktop

desktop = None

# Get SALOME PyQt interface
import SalomePyQt
import libSALOME_Swig

########################################################
# Global variables
########################################################

sgPyQt = SalomePyQt.SalomePyQt()
sg = libSALOME_Swig.SALOMEGUI_Swig()
sgDesktop = sgPyQt.getDesktop()

moduleDesktop   = {}

########################################################
# Internal methods
########################################################

def getStudyId():
    """This method returns the active study ID"""
    return sgPyQt.getStudyId()

def getStudy():
    """This method returns the active study"""

    studyId = _getStudyId()
    study = getStudyManager().GetStudyByID( studyId )
    return study

def getDesktop():
    """This method returns the current TRUST_PLOT2D desktop"""

    global desktop
    return desktop

def setDesktop( studyID ):
    """This method sets and returns TRUST_PLOT2D desktop"""

    global moduleDesktop, desktop

    if not moduleDesktop.has_key( studyID ):
      moduleDesktop[studyID] = CFDesktop( sgPyQt )
      moduleDesktop[studyID].initialize()
    desktop = moduleDesktop[studyID]
    return desktop

################################################
# Callback functions
################################################

def initialize():
    """This method is called when module is initialized. It performs initialization actions"""
    setDesktop( getStudyId() )
    pass

def activeStudyChanged( studyID ):
    """This method is called when active study is changed"""

    setDesktop( getStudyId() )
    pass

def windows():
    """This method is called when module is initialized. It returns a map of popup windows to be used by the module"""
    wm = {}
    return wm

def views():
    """This method is called when module is initialized. It returns a list of 2D/3D views to be used by the module"""
    return []

def activate():
    """This method mimicks SALOME's module activation """
    global desktop
    fv = desktop.showCentralWidget()
    return True

def deactivate():
    """This method is called when module is deactivated"""
    global moduleDesktop, widgetDialogBox
    moduleDesktop[getStudyId()].hideCentralWidget()
