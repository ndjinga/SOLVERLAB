#!/usr/bin/env python
# -*- coding: utf-8 -*-

# %% LICENSE_SALOME_CEA_BEGIN
# see PACKAGESPY/LICENCE file
# %% LICENSE_END


import os
import code
import pprint as PP
import traceback

verbose = False


class InterpreterPython(code.InteractiveInterpreter):
  """
  see https://docs.python.org/2/library/code.html
  """

  def __init__(self, aDict=None):
    """
    InteractiveInterpreter class is not a new-style class
    http://stackoverflow.com/questions/9698614/super-raises-typeerror-must-be-type-not-classobj-for-new-style-class
    """
    if aDict == None:
      aLocalDict = {}
    else:
      aLocalDict = dict(aDict)
    code.InteractiveInterpreter.__init__(self, aLocalDict)
    self._aDict = aLocalDict
    self._stderr = ""
    self._source = ""
    # because uranie ROOT is ModuleFacade
    self._excludeClasses = ["ModuleFacade", "module",
                            "_Feature", "function",
                            "type"]     # TODO what else?

  def compile(self, source):
    """returns compiled source"""
    return code.compile_command(source)

  def runcode_try(self, source):
    """
    try except initial method runcode as python3.9/code.py runcode
    """
    try:
      exec(source, self.locals)
      self._stderr = ""
      # print("OOOOOKKKKK")
    except SystemExit:
      raise
    except Exception as e:
      # print("KKKKKOOOOO %s" % e)
      trace = traceback.format_exc()
      self._stderr = trace
    return

  def runfile(self, aFilePy):
    """set current dir in dir file"""
    with open(aFilePy, "r") as f: code = f.read()
    previousDir = os.getcwd()
    if True: #try:
      newDir = os.path.dirname(os.path.realpath(aFilePy))
      os.chdir(newDir)
      self.runcode(code)
    else: #except:
      pass
    os.chdir(previousDir)

  def write(self, data):
    """
    Write a string to the standard error stream (sys.stderr).
    this derived classes override this to provide the appropriate output handling as needed
    """
    if verbose: print("InterpreterPython.write in stderr: '%s'" % data)
    self._stderr += data

  def getStdErr(self):
    """returns list of data from all write(data) calls as stderr outputs"""
    return self._stderr

  def clearStdErr(self):
    """clear StdErr ouputs"""
    self._stderr = ""

  def getVars(self, excludeClasses=None):
    """get vars in dict of InteractiveInterpreter"""
    #for key, value in d.items():
    if excludeClasses == None:
      exclude = self._excludeClasses
    else:
      exclude = excludeClasses
    newVars = {k: v  for k, v in list(self._aDict.items()) \
               if ( k[0] not in [ "_" ] and \
                    k != "X__source__" and \
                    v.__class__.__name__ not in exclude )  \
               }
    return newVars

  def get__doc__(self):
    try:
      res = self._aDict["__doc__"]
    except:
      res = ""
    return res

  def isOk(self):
    return self.getStdErr() == ""

  def getStrResume(self, source=None):
    res =  "########### result InterpreterPython:\n"
    if source != None:
      #indent source
      indentedSource = "\n".join(["  "+line for line in source.split("\n")])
      res += "\n### code source:\n%s\n" % indentedSource
    err = self.getStdErr()
    if err != '':
      res += "### stderr:\n'%s'\n" % self.getStdErr()
    res += "\n### variables evaluation in code:\n"
    ipvars = self.getVars()
    for k in sorted(ipvars.keys()):
      res += "\n  %s = %s" % (k, PP.pformat(ipvars[k]))
    return res
