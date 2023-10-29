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
import xyzpy.loggingXyz as LOG
import pprint as PP
import settingspy.settings as SETT

logger = LOG.getLogger()

verbose = False

###################################################################
class SolverlabSettings(SETT.Settings):
  """
  may be future link to Qt QSettings
  for future do not forget window environment 
  variables are NOT case sensitive
  user have to write environment variables as syntax ${...}
  
  policy: 
    names setting variables beginning with "_" contains 
    environment variables reference
  """
  
  def __init__(self):
    super(SolverlabSettings, self).__init__()
    self.name = "SolverlabSettings"


#singleton setting for solverlab
_Settings = SolverlabSettings()

def getSettings():
  """used as singleton"""
  return _Settings

def checkEnvVar(val):
  """used as singleton"""
  return _Settings.checkEnvVar(val)

def checkAll():
  """used as singleton"""
  return _Settings.checkAll()

def getVar(name):
  """used as singleton"""
  return _Settings.getVar(name)

def getExpandedVar(name):
  """used as singleton"""
  return _Settings.getExpandedVar(name)

def setEnvVarByDefault(envVar, valueDefault):
  return SETT.setEnvVarByDefault(envVar, valueDefault)
 
def setEnvVar(envVar, value):
  return SETT.setEnvVar(envVar, value)
 
def getSolverlabSysMacrosDir():
  rootdir = getExpandedVar("_SOLVERLABGUI_ROOT_DIR")
  if rootdir is None:
    # raise Exception("user have to set env var 'SOLVERLABGUI_ROOT_DIR'")
    logger.critical("user have to set env var 'SOLVERLABGUI_ROOT_DIR'")
    return None
  aDir = os.path.join(rootdir, "macros")
  if not os.path.isdir(aDir):
    logger.debug("inexisting ${SOLVERLABGUI_ROOT_DIR}/macros dir:\n  '%s'" % aDir)
    #devel uranie for gaudier, current dir
    aDir = os.getcwd()
  logger.debug("getSolverlabSysMacrosDir %s" % aDir)
  return aDir

# policy: 
#   names var begins with "_" contains environment variable var reference
#   theses setting are evaluated LATER,
#   when using USET.getExpandedVar("_SOLVERLABCODE_ROOT_DIR"), for example

#_Settings.setVar("_URANIESYS", "${URANIESYS}")
_Settings.setVar("_SOLVERLABCODE_ROOT_DIR", "${SOLVERLABCODE_ROOT_DIR}")
# _Settings.setVar("_ROOTSYS", "${ROOTSYS}")
#_Settings.setVar("_URANIEWORKDIR", "${URANIEWORKDIR}")
_Settings.setVar("_SOLVERLABGUI_WORKDIR", "${SOLVERLABGUI_WORKDIR}")
#_Settings.setVar("_URANIEGUIDIR", "${URANIEGUIDIR}")
_Settings.setVar("_SOLVERLABGUI_ROOT_DIR", "${SOLVERLABGUI_ROOT_DIR}")


if platform.system() == "Windows":
  _Settings.setVar("editor", "notepad")
  _Settings.setVar("webbowser", "chrome.exe") #r'"C:\Program Files (x86)\Google\Chrome\Application\chrome.exe"')
  _Settings.setVar("terminal", "cmd.exe") #centOs7 and Ubuntu
  _Settings.setVar("difftool", "????") #tkdiff
  _Settings.setVar("gitviewer", "???.") #qgit
  _Settings.setVar("withgit", "no") #qgit
else:
  _Settings.setVar("editor", "pluma")
  _Settings.setVar("webbrowser", "firefox")
  _Settings.setVar("terminal", "gnome-terminal") #centOs7 and Ubuntu
  _Settings.setVar("difftool", "meld") #tkdiff
  _Settings.setVar("gitviewer", "gitk") #qgit
  _Settings.setVar("withgit", "yes") #qgit



#default setting
#parent directory where create etudes directories is current launch directory
#setEnvVarByDefault("${SOLVERLABGUI_WORKDIR}", os.getcwd())
#parent directory where create etudes directories is $HOME/UranieWorkdir
"""
ok, valueDefault = _Settings.checkEnvVar("${HOME}/URANIEWORKDIR")
if ok:
  SETT.setEnvVarByDefault("URANIEWORKDIR", valueDefault)
else: #if not exist $HOME !
  SETT.setEnvVarByDefault("URANIEWORKDIR", os.getcwd())
"""

if verbose:
  print(getSettings())
  print(PP.pformat(getSettings().checkAll()))

