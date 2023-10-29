#!/usr/bin/env python
#-*- coding:utf-8 -*-

# %% LICENSE_SALOME_CEA_BEGIN
# see PACKAGESPY/LICENCE file
# %% LICENSE_END


"""
utilities for best use ConfigParser

see: https://wiki.python.org/moin/ConfigParserExamples
"""

import os
import sys
import pprint as PP

import xyzpy.stringIO as IOX
import xyzpy.loggingXyz as LOG

logger = LOG.getLogger()

# fixing some Python 2-3 code
try:
  import configparser as CPA # 3
  _newStyle = True
except ImportError:
  import ConfigParser as CPA # 2
  _newStyle = False

"""
The ConfigParser :
in the 2.7 version : https://docs.python.org/2/library/configparser.html, but the ExtendedInterpolation() is not present
in the 3.x version : https://docs.python.org/3.2/library/configparser.html
"""

try:
  _ExtendedInterpolation = CPA.ExtendedInterpolation()
except:
  _ExtendedInterpolation = None

_messageDone = [False] # mutable

################################################################
class UtSafeConfigParser(CPA.SafeConfigParser):
  """
  SafeConfigParser with ExtendedInterpolation,
  and __repr__,
  and readDefaultAndUser to merge default and user config
  """

  def __init__(self, *args, **kwargs):
    """
    The SafeConfigParser class has been renamed to ConfigParser
    in Python 3.2. This alias will be removed in future versions.
    Use ConfigParser directly instead.
    (but python 2 yet... for the moment)
    """
    if _newStyle:
      super(UtSafeConfigParser, self).__init__(*args, **kwargs)
    else:
      CPA.SafeConfigParser.__init__(self, *args, **kwargs)
    try:
      self._interpolation = CPA.ExtendedInterpolation()
    except:
      if _messageDone[0] is False:
        logger.warning("UtSafeConfigParser (python2) without ExtendedInterpolation\nin %s" % CPA.__file__)
        # logger.warning("dir(configparser):\n%s" % PP.pformat(dir(CPA)))
        _messageDone[0] = True

  def __repr__(self):
    aDict = self.toDict()
    res = PP.pformat(aDict)
    if len(res) == 0:
      return "UtSafeConfigParser()"
    else:
      return "UtSafeConfigParser(\n %s\n)" % res[1:-1]

  def copy(self):
    return None #TODO

  def toDictTuple(self, verbose=False):
    res = {}
    for sec in self.sections():
      if verbose: print("section %s" % sec)
      opts = []
      for opt in self.options(sec):
        val = self.get(sec, opt)
        if verbose: print("section %s option % s = %s" % (sec, opt, val))
        opts.append((opt, val))
      res[sec] = opts
    return res

  def toOrderedDict(self, verbose=False):
    import collections
    res = collections.OrderedDict()
    for sec in self.sections():
      if verbose: print("section %s" % sec)
      opts = collections.OrderedDict()
      for opt in self.options(sec):
        val = self.get(sec, opt)
        if verbose: print("section %s option % s = %s" % (sec, opt, val))
        opts[opt] = val
      res[sec] = opts
    return res

  def toDict(self, verbose=False):
    res = {}
    for sec in self.sections():
      if verbose: print("section %s" % sec)
      opts = {}
      for opt in self.options(sec):
        val = self.get(sec, opt)
        if verbose: print("section %s option % s = %s" % (sec, opt, val))
        opts[opt] = val
      res[sec] = opts
    return res

  def toCatchAll(self, verbose=False):
    """
    permits class attribute writings,
    (but raise on accentuation and avoid spaces in section names)

    | cfg = UtSafeConfigParser()
    | cfg.readFromStr('''\
    | [General]
    | reporter = tintin
    | ''')
    | config = cfg.toCatchAll()
    | print(config.General.reporter) # -> "tintin"
    """
    from salomepy.catchAll import CatchAll
    res = CatchAll()
    for sec in self.sections():
      if verbose: print("section %s" % sec)
      opts = CatchAll()
      res.__setattr__(sec.replace(" ","_"), opts)
      for opt in self.options(sec):
        val = self.get(sec, opt)
        if verbose: print("section %s option % s = %s" % (sec, opt, val))
        opts.__setattr__(opt, val)
    return res

  def readFromStr(self, aStr, merge=False):
    """allow merge only explicitly"""
    if not merge:
      if not self.isEmpty():
        raise Exception("UtSafeConfigParser not empty")
    stream = IOX.StringIO(aStr)
    name = "default"
    self.readfp(stream, name)

  def writeToStr(self):
    """allow merge only explicitly"""
    stream = IOX.StringIO()
    self.write(stream)
    return stream.getvalue()

  def readDefaultAndUser(self, aStrDefault, aStrUser):
    """merge user overriding origin defaults"""
    if not self.isEmpty():
      raise Exception("UtSafeConfigParser not empty")
    stream = IOX.StringIO(aStrDefault)
    name = "default"
    self.readfp(stream, name)
    stream = IOX.StringIO(aStrUser)
    name = "user"
    self.readfp(stream, name)

  def isEmpty(self):
    if len(self._sections) == 0 and len(list(self._defaults.keys())) == 0:
      return True
    else:
      return False

################################################################
# useful global methods
################################################################
def getConfigFromStr(aStr):
  """simple create config from contents as string"""
  cfg = UtSafeConfigParser()
  cfg.readFromStr(aStr)
  return cfg

def getConfigFromFile(aFile):
  with open(aFile, "r") as f:
    cfg = getConfigFromStr(f.read())
  return cfg


def getConfigFromDefaultAndUserStr(aStrDefault, aStrUser):
  """simple create config from Default an overrrides from User as strings"""
  cfg = UtSafeConfigParser()
  cfg.readDefaultAndUser(aStrUser, aStrDefault)
  return cfg

def toListInt(aStr):
  """transform '1,2,3,...' to list of integer [1, 2, 3, ...]"""
  return [int(i) for i in aStr.split(",")]



