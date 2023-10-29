#!/usr/bin/env python
#-*- coding:utf-8 -*-

# %% LICENSE_SALOME_CEA_BEGIN
# see PACKAGESPY/LICENCE file
# %% LICENSE_END

"""
This file assume DEBUG functionalities use.
Print iradinaGui debug messages in sys.stderr.
Show pretty print debug representation from instances of iradinaGUI classes

| Warning: supposedly show messages in iradinaGUI development phase, not production
| 
| Usage:
| >> import debug as DBG
| >> DBG.write("aTitle", aVariable)        # not shown in production 
| >> DBG.write("aTitle", aVariable, True)  # unconditionaly shown 
"""

import os
import sys
import traceback
import pprint as PP
import iradinapy.coloringIra as COLS

# fixing some Python 2 code but cause pb on unicode
try:
    from StringIO import StringIO
except ImportError:
    from io import StringIO

_debug = [False] #support push/pop for temporary activate debug outputs

_user = os.getenv('USERNAME')
# wambeke is christian at home
_developers = ["christian", "xxwambeke",] # yourself, if you want ...


def isDeveloper():
    """if you are a developer, sometimes you want verbose traces etc."""
    res = _user in _developers
    return res

def indent(text, amount=2, ch=' '):
    """indent multi lines message"""
    padding = amount * ch
    return ''.join(padding + line for line in text.splitlines(True))

def isTypeConfig(var):
    """To know if var is instance from Config/pyconf"""
    typ = str(type(var))
    # print "isTypeConfig" ,type, dir(var)
    if ".pyconf.Config" in typ: return True
    if ".pyconf.Mapping" in typ: return True
    if ".pyconf.Sequence" in typ: return True
    # print "NOT isTypeConfig %s" % typ
    return False
    
def write(title, var="", force=None, fmt="\n#### DEBUG: %s:\n%s\n", fmtEmpty="\n#### DEBUG: %s"):
    """write sys.stderr a message if _debug[-1]==True or optionaly force=True"""
    if _debug[-1] or force:
      tvar = type(var)
      typ = str(tvar)
      if isTypeConfig(var):
        sys.stderr.write(fmt % (title, indent(COLS.toColor(getStrConfigDbg(var)))))
        return
      if var == "": # empty var no write empty line
        sys.stderr.write(fmtEmpty % title)
        return  
      if 'loggingIra.UnittestStream' in typ:
        sys.stderr.write(fmt % (title, indent(var.getLogs())))
        return
      if tvar is not str: # python3 and tvar is not unicode:
        sys.stderr.write(fmt % (title, indent(PP.pformat(var))))
        return
      sys.stderr.write(fmt % (title, indent(var)))
      return
    return

def tofix(title, var="", force=None):
    """
    write sys.stderr a message if _debug[-1]==True or optionaly force=True
    use this only if no logger accessible for classic logger.warning(message)
    """
    fmt = "\n#### TOFIX: %s:\n%s\n"
    write(title, var, force, fmt)

def push_debug(aBool):
    """set debug outputs activated, or not"""
    _debug.append(aBool)

def pop_debug():
    """restore previous debug outputs status"""
    if len(_debug) > 1:
        return _debug.pop()
    else:
        sys.stderr.write("\nERROR: pop_debug: too much pop.")
        return None

def format_color_exception(msg, limit=None, trace=None):
    """
    Format a stack trace and the exception information.
    as traceback.format_exception(), with color
    with traceback only if _debug or isDeveloper())
    """
    etype, value, tb = sys.exc_info()
    if _debug[-1] or isDeveloper():
      res = "<red>" + msg
      if tb:
          res += "<yellow>\nTraceback (most recent call last):\n"
          res += "".join(traceback.format_tb(tb, limit)) #[:-1])
      res += "\n<red>"
      res += "\n".join(traceback.format_exception_only(etype, value))
      return res + "<reset>"
    else: # less verbose
      res = "<red>" + msg # + "<bright>"
      res += "".join(traceback.format_exception_only(etype, value))
      return res+ "<reset>"
      

###############################################
# utilitaires divers pour debug
###############################################

class OutStream(StringIO):
    """
    utility class for pyconf.Config output iostream
    """
    def close(self):
      """
      because Config.__save__ calls close() stream as file
      keep value before lost as self.value
      """
      self.value = self.getvalue()
      StringIO.close(self)
    
class InStream(StringIO):
    """utility class for pyconf.Config input iostream"""
    pass

def getLocalEnv():
    """get string for environment variables representation"""
    res = ""
    for i in sorted(os.environ):
        res += "%s : %s\n" % (i, os.environ[i])
    return res
