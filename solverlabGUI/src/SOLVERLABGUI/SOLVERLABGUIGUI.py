#!/usr/bin/env python
# -*- coding: utf-8 -*-

# %% LICENSE_SALOME_CEA_BEGIN
# see SOLVERLABGUI/LICENSE file
# %% LICENSE_END

import os
import sys
import traceback
import string
import pprint as PP

from PyQt5.QtGui import *
from PyQt5.QtCore import *
from PyQt5.QtWidgets import *

from SOLVERLABDesktop import SOLVERLABDesktop as Desktop
import SOLVERLAB_utils as MUT # as Module_UTils

import salome
# Get SALOME PyQt interface
import SalomePyQt
import libSALOME_Swig

########################################################
# Global variables
########################################################

sgPyQt = SalomePyQt.SalomePyQt()
sg = libSALOME_Swig.SALOMEGUI_Swig()
sgDesktop = sgPyQt.getDesktop()

_NAME = "SOLVERLAB"
# one study
moduleStudy = None
# one desktop
moduleDesktops = {}
# associate ID and action
dict_command = {}

def verbose():
  return False

def log(mess):
  if verbose: MUT.log_print(_NAME+ "GUI" + mess)

########################################################
# Internal methods
########################################################

def getStudy():
    """This method returns the active study"""
    global moduleStudy, moduleDesktops
    current = str(salome.myStudy)
    if moduleStudy is None:
       moduleStudy = current  # first study detected copy objref_study as str
       # first study as '<SALOMEDS._objref_Study object at 0x59398a8>' for example
       log(" getStudy first activation on study %s" % current)
    else:
       if moduleStudy != current: # user have done menu "new study"
          log(" WARNING: problematic new study detected: Relaunch Salome is the best recommended practice\n")
          # problematic as wrapped C/C++ object of type QMainWindowForLog has been deleted and etc on menu and actions/slots
          moduleDesktops = {}
          setDesktop()
          moduleStudy = current
    # log(" getStudy %s %s" % (salome.myStudy, moduleStudy))
    return salome.myStudy

def getDesktop():
    """This method returns the current GUI desktop"""
    MUT.log_print(" getDesktop")
    return sgPyQt.getDesktop()

def setDesktop():
    """This method sets and returns GUI desktop"""
    global moduleDesktops, dict_command
    log(" setDesktop()")
    study = str(salome.myStudy)
    if study not in moduleDesktops:
        log(" create Desktop")
        desk = Desktop( sgPyQt, sg )
        moduleDesktops[study] = desk
        desk.set_dict_command( dict_command )
        # log(" set dict_command", dict_command)
    return moduleDesktops[study]

def resetDesktop():
    """This method sets and returns GUI desktop widgets (when there is new study,for example)"""
    log(" TODO resetDesktop()")

def incObjToMap(m, id):
    """This method increments the object counter in the map"""
    if id not in m:
        m[id] = 0
    m[id] += 1

def getSelection():
    """This method analyses selection"""
    selcount = sg.SelectedCount()
    seltypes = {}
    for i in range(selcount):
        incObjToMap( seltypes, getObjectID(getStudy(), sg.getSelected(i)) )
    return selcount, seltypes

################################################
# Callback functions
################################################

def initialize():
    """This method is called when module is initialized. It performs initialization actions"""
    #for i in dir(SalomePyQt): log("SalomePyQt",i
    log(" initialize()")
    setDesktop()

def windows():
    """This method is called when module is initialized. It returns a map of popup windows to be used by the module"""
    #SalomePyQt WT_LogWindow,  WT_ObjectBrowser, WT_PyConsole, WT_User
    log(" windows()")
    wm = {}
    wm[SalomePyQt.WT_ObjectBrowser] = Qt.LeftDockWidgetArea
    # wm[SalomePyQt.WT_PyConsole] = Qt.BottomDockWidgetArea
    return wm

def views():
    """
    for centralwidget, this method is called when module is initialized.
    It returns a list of 2D/3D views to be used by the module
    """
    log(" views()")
    return []

def createPreferences():
    """This method is called when module is initialized. It exports module preferences"""
    pass

def activate():
    """This method is called when module is initialized. It returns True if activating is successfull, False otherwise"""
    global moduleDesktops
    log(" activate()")
    study = str(getStudy())
    if study in moduleDesktops:  # may be change with new study
      moduleDesktops[study].activateDesktop()
    else:
      log(" activate(), no existing desktop")
    return True

def viewTryClose(widID):
    """
    if ViewClosable, delete object,
    generate RuntimeError: wrapped C/C++ object of type QMainWindowForLog has been deleted
    so not ViewClosable, hide it
    """
    log(" viewTryClose()")
    sgPyQt.setViewClosable(widID, False)
    sgPyQt.setViewVisible(widID, False)

def deactivate():
    """This method is called when module is deactivated"""
    global moduleDesktops
    log(" deactivate()")
    study = str(getStudy())
    if study in moduleDesktops:  # may be change with new study
      moduleDesktops[study].deactivateDesktop()
    else:
      print(" deactivate(), no existing desktop")

def activeStudyChanged():
    """This method is called when active study is changed, and when salome desktop is activated"""
    log(" activeStudyChanged %s" % getStudy())

def createPopupMenu(popup, context):
    """
    for salome object browser and centralwidget views,
    this method is called when popup menu is invocked
    """
    log(" createPopupMenu()")

def OnGUIEvent(commandID):
    """This method is called when a GUI action is activated"""
    global dict_command
    log(" OnGUIEvent(): command id %s" % commandID)
    if commandID in dict_command:
       dict_command[commandID]()
    else:
       log("The command id %s is not implemented" % commandID)

def preferenceChanged(section, setting):
    """This method is called when module's preferences are changed"""
    log(" referenceChanged(): section %s" % section)
    pass

def activeViewChanged(viewID):
    """This method is called when active view is changed"""
    log(" activeViewChanged(): id %s" % viewID)
    pass

def viewCloned(viewID):
    """This method is called when active view is cloned"""
    log(" viewCloned(): id %s" % viewID)
    pass

def viewClosed(viewID):
    """This method is called when active view viewClosed"""
    log(" viewClosed(): id %s" % viewID)
    pass

def engineIOR():
    """This method is called when study is opened. It returns engine IOR"""
    res = getEngineIOR()
    log(" engineIOR(): %s" % res)
    return res
