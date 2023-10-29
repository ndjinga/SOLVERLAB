#!/usr/bin/env python
# -*- coding: utf-8 -*-

# %% LICENSE_SALOME_CEA_BEGIN
# see PACKAGESPY/LICENCE file
# %% LICENSE_END


"""
assume singleton of a global variable set from VariablesInterpreterPythonXyz
as copy of real value
useful for all caller
"""

import os
import pprint as PP

import xyzpy.loggingXyz as LOG
logger = LOG.getLogger()

_GlobalVariables = {} # empty dict for named variables
verbose = False

def getGlobalVariables(treePyName=".GlobalVariables"):
  try:
    return _GlobalVariables[treePyName]
  except:
    if verbose:
      logger.warning("tag '%s' not found" % treePyName)
    return None

def setGlobalVariables(treePyName, value):
  _GlobalVariables[treePyName] = value
  if verbose:
    logger.warning("setGlobalVariables %s\n%s" % (treePyName, PP.pformat(_GlobalVariables)))



