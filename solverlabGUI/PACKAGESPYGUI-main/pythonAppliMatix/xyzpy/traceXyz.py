#!/usr/bin/env python
# -*- coding: utf-8 -*-

# %% LICENSE_SALOME_CEA_BEGIN
# see PACKAGESPY/LICENCE file
# %% LICENSE_END

"""
Easy meaning light write operations of log/trace and error and warning and debug
"""

import string
import time
import os
import sys
import subprocess

verbose = False

class TraceXyz(object):
  """
  EZ meaning light class managing write operations of log/trace and error and warning and debug
  used for display foreground and background trace of executing subprocess
  in one or more standalone (...or centralwidget) salome window as qTextEditForLog instances
  with color and time display

  :usage:
  
  >>> from xyzpy.traceXyz import TraceXyz
  >>> CT = TraceXyz(prefix = "combs", windowTitle="job combs of the day")
  >>> ...
  >>> CT.debug("a simple debug message")
  >>> CT.trace("a simple trace message")
  >>> CT.warning("a simple warning message")
  >>> CT.error("a important error message")
  
  """

  def __init__( self, prefix="-", windowTitle="Process Trace"):
    self._controller = None #to get (or not) a logwindow for cout traces...
    self._Trace = True
    self._Error = True
    self._Warning = True
    self._Debug = False
    self.prefix = prefix #no prefix name by default
    self.windowTitle = windowTitle
    self._logWindow = None
    self._lastTrace = "" #to avoid write "process xxx %" if same time
    self._forcePrint = False #to force logCmdWidget print PLUS console print
    #self.onConsole = onConsole #on window plus unconditionnaly on console
    self._noGui = False
  
  #def setOnConsole(onConsole=True):
  #  self.onConsole = onConsole
    
  def strList(self, aList):
    """format %.4e for floats in list"""
    if type(aList) == list:
      res = ""
      for f in aList:
        if type(f) == float:
          res += "%.4e, " % f
        else:
          res += "%s, " % f
      return "[%s]" % res[:-2]
    else: #default classical way
      return str(aList)

  def setNoGui(self):
    self._noGui = True

  def setGui(self):
    self._noGui = False

  def setController(self, controller):
    self._controller = controller
    if verbose: print("TraceXyz.setController", self._controller)

  def coutStd( self, mess ):
    if self._controller == None or self.forcePrint == True:
      print(mess) #simple stdout print
      sys.stdout.flush()
    else:
      logCmdWidget = self._controller.getLogCMDWidget()
      if logCmdWidget != None:
        logCmdWidget.insertText( mess )
      else:
        print(mess) #simple stdout print

  def coutInStandalone( self, mess ):
    if self._noGui:
      print(mess) #simple stdout print
      sys.stdout.flush()
      return

    if self._logWindow == None:
      try:
        import salomepy.qTextEditForLog as QTL
        self._logWindow =  QTL.QTextEditForLog()
        self._logWindow.setWindowTitle(self.prefix)
      except:
        self.setNoGui() #no way to GUI...
        self.coutInStandalone(mess)
        return
        
    self._logWindow.show()
    color = "Black"
    if "ERROR:" in mess: color = "Red"
    if "WARNING:" in mess: color = "Blue"
    if "DEBUG:" in mess: color = "Magenta"
    self._logWindow.insertTextColor( mess, color)

  def cout( self, mess ):
    self.coutInStandalone( mess )
    if self._forcePrint:
      print(mess)

  def trace( self, txt ):
    """
    trace computation
    filter ': process xx %' of same time
    """
    if not self._Trace: return #no trace

    mess = "%s: %s %s" % (self.prefix, time.strftime('%X:'), txt)
    aStr = ": process "
    if aStr in mess:
      if aStr in self._lastTrace:
        aSpl = mess.split(aStr)[0]
        if aSpl == self._lastTrace.split(aStr)[0]:
          if " 100 %" not in mess:
            return #no print process xx% of same time, but 100 %
    self._lastTrace = mess
    self.cout(mess)

  def info( self, txt, withStack=False ):
    """idem CT.trace"""
    self.trace( txt, withStack=withStack )

  def error( self, txt ):
      """error computation"""
      if self._Error:
          mess = "ERROR: %s: %s %s" % (self.prefix, time.strftime('%d/%m/%y %X:'), txt)
          self.cout(mess)

  def warning( self, txt ):
      """warning computation"""
      if self._Warning:
          mess = "WARNING: %s: %s %s" % (self.prefix, time.strftime('%d/%m/%y %X:'), txt)
          self.cout(mess)

  def debug( self, txt ):
      """debug computation"""
      if self._Debug:
          mess = "DEBUG: %s: %s %s" % (self.prefix, time.strftime('%d/%m/%y %X:'), txt)
          self.cout(mess)

"""
here create one instance, use for one common widget,
example for only one compute job at a time in foreground
attach others instance TraceXyz() to others
compute combs jobs in background to avoid mixing traces
"""
GeneralTrace = TraceXyz(windowTitle="General Trace")

