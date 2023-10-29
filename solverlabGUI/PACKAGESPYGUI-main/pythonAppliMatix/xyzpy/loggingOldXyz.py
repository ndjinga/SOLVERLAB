#!/usr/bin/env python
# -*- coding: utf-8 -*-

# %% LICENSE_SALOME_CEA_BEGIN
# see PACKAGESPY/LICENCE file
# %% LICENSE_END

"""
Logger for packagespy using logging package
OLD instantiation - before MATIX V3

Manage write operations of log/trace of Xyz package:

- used for display foreground trace of executing process
- using logging python package
- outputs in console or centralwidget salome window (as qTextEditForLog) or other...
- outputs trace of file.py name and number of line of call ...

:usage:

>>> import xyzpy.loggingXyz as LOG
>>> logger = LOG.getLogger()
>>> ...
>>> logger.debug("a simple debug message")
>>> logger.info("a simple trace message")
>>> logger.warning("a simple warning message")
>>> logger.error("a important error message")

"""

import os
import logging as LOGI
import pprint as PP

print("loggingOldXyz Import - OLD instantiation")

def verbose():
  return False

_STEP = LOGI.INFO - 1 # step level is just below INFO
_TRACE = LOGI.INFO - 2 # trace level is just below STEP

LOGI.addLevelName(_STEP, "STEP")
LOGI.addLevelName(_TRACE, "TRACE")
LOGI.STEP = _STEP
LOGI.TRACE = _TRACE

if verbose():
  print('dir(LOGI)\n%s' % PP.pformat(dir(LOGI)))

if False:
  # skip initialisations as logger set yet somewhere else
  _logger = LOGI.getLogger("SomewhereElseDefaultLogger")
  _setLevelDone = [True]
  _setconsoleHandlerDone = [True]
  _logger.info("------logger XyzLogger set as XyzDefaultLogger------")
else:
  # initialisations as logger not set yet somewhere else
  _logger = LOGI.getLogger("XyzLogger")
  _setLevelDone = [False]
  _setconsoleHandlerDone = [False]

__pushpopLevels__ = []

if verbose():
  print('dir(_logger)\n%s' % PP.pformat(dir(_logger)))
  print("logging.__file__", LOGI.__file__)
  print("logginXyz _logger.name", _logger.name)


levels = {
  "CRITICAL": LOGI.CRITICAL,
  "ERROR": LOGI.ERROR,
  "WARNING": LOGI.WARNING,
  "INFO": LOGI.INFO,
  "STEP": _STEP, # step level is just below INFO
  "TRACE": _TRACE, # trace level is just below STEP,
  "DEBUG": LOGI.DEBUG
}

_levelNames = {
  LOGI.CRITICAL: "CRITICAL",
  LOGI.ERROR: "ERROR",
  LOGI.WARNING: "WARNING",
  LOGI.INFO: "INFO",
  LOGI.STEP: "STEP",
  LOGI.TRACE: "TRACE",
  LOGI.DEBUG: "DEBUG"
}

_knownLevels = "CRITICAL ERROR WARNING INFO STEP TRACE DEBUG".split()
_knownLevelsStr = "[%s]" % "|".join(_knownLevels)

_MISSED_GETENV = []

def close():
  print("close TODO")
  return

_logger.trace = _logger.info # more simple for the moment
_logger.step = _logger.info # more simple for the moment
_logger.close = close


_logger.setLevel(LOGI.INFO)
#_logger.info('initial info message')

if verbose(): print("dir(logging):\n", dir(LOGI))

#################################################################
def indent(msg, nb=11, car=" "):
  """indent nb car (spaces) multi lines message except first one"""
  s = msg.split("\n")
  res = ("\n" + car * nb).join(s)
  return res

#################################################################
def indentUnittest(msg, prefix="            "): # older version prefix=" | "
  """
  indent multi lines message except first one with prefix.
  prefix default is designed for less spaces for size logs files
  and keep logs human eye readable
  """
  s = msg.split("\n")
  res = ("\n" + prefix).join(s)
  return res

########################################################################

def filterLevel(aLevel):
  """
  filter levels logging values from firsts characters levels.
  No case sensitive.

  | example:
  | 'i' -> 'INFO'
  | 'cRiT' -> 'CRITICAL'
  """
  aLev = aLevel.upper()
  knownLevels = _knownLevels
  maxLen = max([len(i) for i in knownLevels])
  for i in range(maxLen):
    for lev in knownLevels:
      if aLev == lev[:i]:
        # DBG.write("filterLevel", "%s -> %s" % (aLevel, lev))
        return lev
  msg = "Unknown level '%s', accepted are: %s" % (aLev, ", ".join(knownLevels))
  raise Exception(msg)

#################################################################
def tryLogger():
  _logger.debug('try debug message')
  _logger.info('try info message')
  _logger.warning('try warning message')
  _logger.error('try error message')
  _logger.critical('try critical message')

#################################################################
def setLevel(level=None):
  if verbose():
    print("LOG.setLevel", level)
  _setLevelDone[0] = True
  if level == None:
    consoleHandler.setLevel(LOGI.WARNING)  # default, no message
    return
  try:
    consoleHandler.setLevel(int(level))
    _logger.setLevel(int(level))
    return
  except:
    pass
  if level in list(levels.keys()):
    consoleHandler.setLevel(levels[level])
    _logger.setLevel(levels[level])
    print("bonjour %s" % levels[level])
    return
  if 'Level ' in level:
    try:
      consoleHandler.setLevel(int(level[6:]))
      _logger.setLevel(int(level[6:]))
      return
    except:
      pass

  consoleHandler.setLevel(LOGI.INFO)  # default
  _logger.setLevel(LOGI.INFO)
  _logger.warning(("Bad logger level: '%s' try %s"), level, list(levels.keys()))
  return

#################################################################
def getLogger():
  if _setLevelDone[0] == False: setLevel("WARNING")
  return _logger

#################################################################
def getLoggerLevelName():
  try:
    return LOGI.getLevelName(consoleHandler.level)
    # return LOGI._levelNames[consoleHandler.level]
    # return LOGI._levelNames[logger.getEffectiveLevel()]
  except:
    try:
      return _levelNames[consoleHandler.level]
    except:
      return str(consoleHandler.level)

#################################################################
def pushLevel(level):
  __pushpopLevels__.append(getLoggerLevelName())
  setLevel(level)

#################################################################
def popLevel(warning=True):
  if len(__pushpopLevels__) > 0:
    level = __pushpopLevels__.pop()
    setLevel(level)
  else:
    if warning: _logger.warning("Pop logger level empty")

#################################################################
def getenv(name):
  """
  signal only one time missing env var
  returns empty string in this case, not None
  """
  res = os.getenv(name)
  if res == None:
    if name not in _MISSED_GETENV:  # warning problem only one time
      _MISSED_GETENV.append(name)
      _logger.warning("missed environment variable '%s'" % name)
    return ""  # avoid raise error, problems in os.path.join() AttributeError: 'NoneType' object has no attribute 'endswith'
  return res

#################################################################
class UnittestFormatter(LOGI.Formatter, object): # object force new-style classes in logging 0.5.0.5 python 2.6
  """
  this formatter prefixes level name and indents all multi lines messages
  """
  def format(self, record):
    res = super(UnittestFormatter, self).format(record)
    # res = indentUnittest(res, prefix="           ")
    res = indentUnittest(res)
    return res

#################################################################
class DefaultFormatter(LOGI.Formatter, object): # object force new-style classes in logging 0.5.0.5 python 2.6
  """
  this formatter prefixes level name and indents all multi lines messages
  """
  def format(self, record):
    res = super(DefaultFormatter, self).format(record)
    # res = indent(res,)
    res = indentUnittest(res)
    return res

#################################################################
class DefaultFormatter_as_sat(LOGI.Formatter, object): # object force new-style classes in logging 0.5.0.5 python 2.6
  """
  this formatter prefixes level name and indents all messages but INFO stay "as it" (as SAT does)
  """
  def format(self, record):
    # print "", record.levelname #type(record), dir(record)
    # nb = len("2018-03-17 12:15:41 :: INFO     :: ")
    if record.levelname == "INFO":
      res = record.getMessage()
    else:
      res = super(DefaultFormatter, self).format(record)
      # res = indent(res,)
      res = indentUnittest(res, prefix=" | ")
    return res


#if len(LOGI._handlerList)==0:
if _setconsoleHandlerDone[0] == False:
  if verbose(): print("WARNING: add my StreamHandler to logging")

  # create console handler with my log level
  consoleHandler = LOGI.StreamHandler()

  # create formatter and add it to the handlers
  user = os.getenv("USERNAME")
  if user in ["xxcvw", "xxchristian", "xxwambeke", "ym268439"]:
    # formatter = UnittestFormatter('%(levelname)-8s : %(name)-15s : %(filename)-15s (%(lineno)3d) : %(funcName)-10s : %(message)s')
    formatter = UnittestFormatter('%(levelname)-8s : %(filename)-15s (%(lineno)3d) : %(message)s')
    consoleHandler.setLevel(LOGI.DEBUG)
  else:
    consoleHandler.setLevel(LOGI.WARNING)
    formatter = DefaultFormatter('%(levelname)-8s : %(message)s')
  consoleHandler.setFormatter(formatter)
  # add the handlers to the logger
  _logger.addHandler(consoleHandler)
  _setconsoleHandlerDone[0] = True

  if True: #TODO verbose():
    _logger.info("set new consoleHandler in logger '%s'" % _logger.name)
    tryLogger()

  # this creates big trouble in salome GUI...
  # from salomepy.onceQApplication import OnceQApplication
  #_logger.warning('logger create a QApplication')
  # app = OnceQApplication([])
