#!/usr/bin/env python
#-*- coding:utf-8 -*-

#  Copyright (C) 2010-2023  CEA
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

"""
This file is the main API for config configparser for solverlabGUI
"""

import os
import platform
import shutil
import fnmatch
import sys

import debogpy.debug as DBG
import returncodepy.returnCode as RCO
import solverlabpy.utilsSvl as UTS
import solverlabpy.dateTime as DATT

import configparserpy.configParserUtils as CPAU

_verbose = False

_configUserStr = """\
# User config file solverlabGUI_user.cfg
# user may modify here configuration at his convenience...
# original configuration values are in config file solverlabGUI_default.cfg
# see: https://wiki.python.org/moin/ConfigParserExamples

# ... at your own risks ...
# WARNING: names of values are in lowercase evaluation

# example...

[MainWindow]
title = Solverlab GUI
sizex = 900
sizey = 600

"""

_configDefaultStr = """\
# User config file solverlabGUI_default.cfg
# user find here default configuration values
# see: https://wiki.python.org/moin/ConfigParserExamples

# ALL USER MODIFICATIONS HERE WILL BE REMOVED
# WARNING: names of values are in lowercase evaluation

[General]
# usermode implies what data to show: advanced, simple (for now)
usermode = simple

[MainWindow]
title = Solverlab GUI
sizex = 900
sizey = 600
color_treeview_base = 190, 190, 205
color_treeview_text = 0, 0, 0

"""

"""
#TODO future

[DefaultEnvironment]
# TODO not used yet
# solverlab code root dir (.../solverlab)
solverlab_root_dir = default
# solverlab GUI working dir (.../SOLVERLABGUI_WORKDIR)
solverlabgui_workdir = default
# solverlab GUI LOGS dir  (.../SOLVERLABGUI_WORKDIR/LOGS)
solverlabgui_logdir = default

[Logger]
logdir = LOG

"""

_mainConfig = None # singleton for a general config
_mainConfigManager = None # singleton for a general config manager
_SOLVERLABGUI_WORKDIR = None

# import methods for EZ direct use from configSvl
toListInt = CPAU.toListInt

class ConfigManager(object):
  """
  Manages the read/write of config files of solverlabGUI,
  and merges if useful.

  | file solverlabGUI_user.cfg
  | file solverlabGUI_default.cfg
  """

  def __init__(self, runner=None):

    if runner is not None:
      self.runner = runner
      # self.logger = runner.getLogger()
      self.options = runner.getOptions()
    else:
      self.runner = None
      # self.logger = None
      self.options = None

    self.fileUser = "solverlabGUI_user.cfg"
    self.fileDefault = "solverlabGUI_default.cfg"
    # better store string as ready to create new config instance(s)
    # and readable for debug
    self.configUserStr = None # done one time
    self.configDefaultStr = None # done one time
    self._configUserStr = _configUserStr
    self._configDefaultStr = _configDefaultStr

  def getWorkdir(self):
    global _SOLVERLABGUI_WORKDIR
    if self.options is not None:
      return self.options.workdir
    else:
      if _SOLVERLABGUI_WORKDIR is None:
        default = "${HOME}/SOLVERLABGUI_WORKDIR"
        res = os.getenv("SOLVERLABGUI_WORKDIR", default)
        _SOLVERLABGUI_WORKDIR = os.path.realpath(os.path.expandvars(res))
        if _verbose: print("configSvl default workdir %s" % _SOLVERLAB_WORKDIR)
      return _SOLVERLABGUI_WORKDIR

  def getRealPath(self, name):
    res = os.path.join(self.getWorkdir(), name)
    return os.path.realpath(res)

  def checkFileExist(self, filename):
    """filename as name relative to workdir"""
    return os.path.isfile(self.getRealPath(filename))

  def assertUserDefaultFiles(self):
    """
    if inexisting, create config files user and default.
    relative to workdir
    """
    aFile = self.getRealPath(self.fileUser)
    if not os.path.isfile(aFile):
      # inexisting, create it
      open(aFile, "w").write(self._configUserStr)

    # unconditionnaly override it
    aFile = self.getRealPath(self.fileDefault)
    with open(aFile, "w") as f:
      f.write(self._configDefaultStr)

  def _getConfig(self, fileName):
    """create config from a file name, reading file, new instance config"""
    cfg = CPAU.getConfigFromStr(open(aFile).read())
    return cfg

  def getUserConfig(self):
    """new instance user config"""
    if self.configUserStr == None:
      aFile = self.getRealPath(self.fileUser)
      self.configUserStr = open(aFile, "r").read()
    cfg = CPAU.getConfigFromStr(self.configUserStr)
    return RCO.ReturnCode("OK", "get user config", cfg)

  def getDefaultConfig(self):
    """new instance default config"""
    if self.configDefaultStr == None:
      aFile = self.getRealPath(fileDefault)
      self.configDefaultStr = open(aFile, "r").read()
    cfg = CPAU.getConfigFromStr(self.configDefaultStr)
    return RCO.ReturnCode("OK", "get default config", cfg)

  def getMainConfig(self):
    """
    main as merged config of default plus overrides of user
    new instance main config
    """
    self.assertUserDefaultFiles()
    aFile = self.getRealPath(self.fileUser)
    self.configDefaultStr = open(aFile, "r").read()
    aFile = self.getRealPath(self.fileDefault)
    self.configUserStr = open(aFile, "r").read()
    cfg = CPAU.getConfigFromDefaultAndUserStr(self.configDefaultStr, self.configUserStr)
    DBG.write("getMainConfig", cfg.toDict())
    return RCO.ReturnCode("OK", "get main config", cfg)

  def setMainConfig(self, cfg):
    """set a user main config as global, is user choice"""
    _mainConfig[0] = cfg


#########################################################################
# user mode (advanced or else) to show or hide some data in treeview
#########################################################################

# list of hidden items by association with value of item.getTreePyName()

_advanced = [] # [] as nothing hidden

_simple = """
Controls
dataInformations.release
# Components*.Density
# Materials*.TargetConcentration
# Materials*.IsVacuum
""".split()
# fnmatch pattern [seq] matches any character in seq

_modesNames = "simple advanced".split()

_modes = {
  "simple": _simple,
  "advanced": _advanced,
}

_currentMode = ["simple"] # list as global mutable
if _verbose: print("set _currentMode", _currentMode[0])

"""
def getUserName():
  res = os.getenv('USERNAME')
  if res == None:
    res = os.getenv('USER')
  if res is None:
    raise Exception("can't get user name in env var 'USER' or 'USERNAME'")
  return res

# set as default if not wambeke or developper
if getUserName() in "wambeke christian".split():
  _currentMode = ["advanced"] # list as global mutable
else:
  _currentMode = ["simple"] # list as global mutable
"""

def getExistingModes():
  # DBG.write("congigSvl.isHidden modes", _modes, DBG.isDeveloper())
  return _modesNames

def setCurrentMode(modeName):
  if not modeName in getExistingModes():
    raise Exception("unknown mode '%s'" % modeName)
  _currentMode[0] = modeName

def getCurrentMode():
  res = _currentMode[0]
  if _verbose: print("getCurrentMode", res)
  return res

def isHidden(item, nameAttr=None, modeName=None):
  """
  avoid Components[*].Density because fnmatch pattern [seq] matches any character in seq.
  use Components*.Density instead
  """
  if modeName is None:
    mode = _currentMode[0]
  else:
    mode = modeName
  name = item.getTreePyName()
  if nameAttr is not None:
    name += "." + nameAttr
  hidden = _modes[mode]
  res = False
  for i in hidden:
    if fnmatch.fnmatch(name, "*" + i):
      # if _verbose: print("fnmatch hidden ", name, "*" + i)
      res = True
      break
    # else:
    #  if "Density" in name:print("not fnmatch hidden ", name, "*" + i)
  # DBG.write("congigSvl.isHidden '%s' %s" % (name, res), "")
  return res


def getMainConfig():
  global _mainConfig
  global _mainConfigManager
  if _mainConfig is None:
    _mainConfigManager = ConfigManager()
    rc = _mainConfigManager.getMainConfig() # current config from options.workdir/xxx.cfg
    if not rc.isOk():
      raise Exception("Problem getMainConfig : %s" % rc)
    _mainConfig = rc.getValue()
  if _verbose: print("configSvl.getMainConfig()\n%s" % _mainConfig)
  return _mainConfig  # set as main singleton ConfigManager.getMainConfig()

def getMainConfigCatchAll():
  res = getMainConfig()
  return res.toCatchAll() # have to be set in ConfigManager.getMainConfig()
