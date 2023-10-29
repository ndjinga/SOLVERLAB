#!/usr/bin/env python
# -*- coding: utf-8 -*-

# %% LICENSE_SALOME_CEA_BEGIN
# see PACKAGESPY/LICENCE file
# %% LICENSE_END


import os
import xyzpy.loggingXyz as LOG
import pprint as PP

logger = LOG.getLogger()

verbose = False

###################################################################
class Settings(object):
  """
  may be future link to Qt QSettings
  for future do not forget window environment variables are NOT case sensitive
  user have to write environment variables as syntax ${...}
  
  policy: names setting variables beginning with "_" contains environment variables reference
  """
  
  def __init__(self):
    self.name = "Settings"
    self._vars = {}

  def setVar(self, name, value, Verbose = True):
    """no control, user choice with checkAll()"""
    if name in list(self._vars.keys()):
      if Verbose: logger.warning("override settings name var: '%s': '%s' -> '%s'" % (name, self._vars[name], value))
    self._vars[name] = value  

  def getVar(self, name, Verbose = True):
    """return NOT expanded value of var name"""
    if name in list(self._vars.keys()):
      return self._vars[name]
    else:
      if Verbose: logger.warning("problem unknown settings name var: '%s'" % name)
      return None 
    
  def getExpandedVar(self, name):
    """
    returns expanded value of var name, 
    or None if inexisting or problem in expand
    """
    if name in list(self._vars.keys()):
      ok, res = self.checkEnvVar(self._vars[name])
      if ok: 
        return res
      else:
        return None
    else:
      logger.warning("problem unknown settings name var: '%s'" % name)
      return None 

  def getRealPath(self, aPathWithEnvVar):
    """
    resolve file path env variable as ${HOME}/toto etc... 
    with os.path.expandvars interpretation of env var
    """
    res, why = self.checkEnvVar(self, aPathWithEnvVar)
    if not ok:
      raise Exception("Problem interpretation env var: %s" % why)
    try:
      res = os.path.expandvars(aPathWithEnvVar)
    except:
      raise Exception("Problem interpretation of path: '%s'" % aPathWithEnvVar)
    
    res = os.path.realpath(res)
    if verbose: logger.info("getRealPath: '%s' -> '%s'" % (aPathWithEnvVar, res))
    return res

  def checkEnvVar(self, aStrWithEnvVar):
    """check all env vars contained as syntax ${...} defined in environ"""
    try:
      res = os.path.expandvars(aStrWithEnvVar)
    except:
      raise Exception("Problem interpretation of env var: '%s'" % aStrWithEnvVar)
    if "$" in res: #unknown env var stay in res
      return (False, '%s -> %s' % (aStrWithEnvVar, res))
    return (True, res)

  def __repr__(self):
    return "%s:\n%s" % (self.name, PP.pformat(self._vars, indent=2))

  def checkAll(self):
    """
    expand var in settings
    return (ok, aDict) ok is False or True
    aDict is settings as aDict[key] = (value, interpretedValue)
    """
    okres = True
    result = {}
    for k, v in list(self._vars.items()):
      if '_' == k[0]:
        ok, why = self.checkEnvVar(v)
        if not ok: 
          okres = False
          why = None
      else:
        ok, why = True, v
      result[k] = (v, why)
    return okres, result

#singleton use
_Settings = Settings()

def getSettings():
  """used as singleton"""
  return _Settings

def getVar(name):
  """used as singleton"""
  return _Settings.getVar(name)

def getExpandedVar(name):
  """used as singleton"""
  return _Settings.getExpandedVar(name)

def setEnvVarByDefault(envVar, valueDefault):
   """
   only for single environ variables:
   envVar='HOME' for ${HOME} or $HOME
   valueDefault have to be expanded (i.e. without '$')
   """
   if "$" in envVar:
     raise Exception("Problem envVar without '$' please : '%s'" % envVar)
   if "$" in valueDefault:
     raise Exception("Problem valueDefault without '$' please : '%s'" % valueDefault)
   dollarEnvVar = "${%s}" % envVar
   ok, res = _Settings.checkEnvVar(dollarEnvVar) #use singleton for ease
   #print 'setEnvVarByDefault', dollarEnvVar, ok, res
   if not ok:
     logger.warning("$%s='%s'" % (envVar, valueDefault))
     os.environ[envVar] = valueDefault
 
def setEnvVar(envVar, value):
   """
   with message warning if change
   """
   logger.debug("setEnvVar %s=%s" % (str(envVar), str(value)))
   if "$" in envVar:
     raise Exception("Problem envVar without '$' please : '%s'" % envVar)
   if "$" in value:
     raise Exception("Problem valueDefault without '$' please : '%s'" % value)
   oldValue = None
   try:
     oldValue = os.environ[envVar]
   except:
     oldValue = None
   if oldValue != value:
     logger.warning("\n  $%s='%s' ->\n  $%s='%s'" % (envVar, oldValue, envVar, value))
     os.environ[envVar] = value
   return
 


