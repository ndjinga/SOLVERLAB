#!/usr/bin/env python
# -*- coding: utf-8 -*-

# %% LICENSE_SALOME_CEA_BEGIN
# see PACKAGESPY/LICENCE file
# %% LICENSE_END


import os
import traceback
from datetime import datetime
import xyzpy.loggingXyz as LOG

logger = LOG.getLogger()

verbose=False

_incrementOnNames = {}

class IncrementOnName(object):
  def __init__(self, anInstance, aPrefix=None):
    """
    for future: self.affectedNames could avoid duplicate existing names (yet in geom for example)
    """
    self.affectedNames = []
    if aPrefix == None:
      self.prefix = str(anInstance.__class__.__name__) + "_"
    else:
      self.prefix = str(aPrefix)
    self.no = 0
  
  def getIncrement(self, aPrefix=None):
    while True:
      if aPrefix == None:
        res = self.prefix + str(self.no)
      else:
        res = aPrefix + str(self.no)
      self.no += 1
      if res not in self.affectedNames: break
    self.affectedNames.append(res)
    return res

def getIncrement(anInstance, aPrefix=None):
  """
  return string like aPrefix_xx
  with xx is an integer with increment on successive calls
  """
  global _incrementOnNames
  className = anInstance.__class__.__name__
  if className not in _incrementOnNames:
    _incrementOnNames[className] = IncrementOnName(anInstance, aPrefix)
  return _incrementOnNames[className].getIncrement(aPrefix)

def resetIncrement(anInstance):
  global _incrementOnNames
  try:
    className = anInstance.__class__.__name__
    del(_incrementOnNames[className])
  except:
    pass

def setAffectedNames(anInstance, affectedNames):
  global _incrementOnNames
  className = anInstance.__class__.__name__
  if className not in _incrementOnNames:
    _incrementOnNames[className] = IncrementOnName(anInstance)
  ic = _incrementOnNames[className]
  for i in affectedNames:
    if i not in ic.affectedNames:
      ic.affectedNames.append(i)
  return

def resetAllIncrements():
  global _incrementOnNames
  _incrementOnNames = {}
