#!/usr/bin/env python
# -*- coding: utf-8 -*-

# %% LICENSE_SALOME_CEA_BEGIN
# see PACKAGESPY/LICENCE file
# %% LICENSE_END


"""
http://stackoverflow.com/questions/874815/how-do-i-get-real-time-information-back-from-a-subprocess-popen-in-python-2-5?lq=1
"""

import threading
from PyQt5 import QtGui, QtWidgets

verbose = False

def toHexAll0x(aStr):
  res = ""
  for x in aStr: res += " " + hex(int(ord(x)))
  return res[1:]

def toHexAll(aStr):
  import binascii
  res = binascii.hexlify(aStr.encode())
  return str(res)[2:-1] # from b'xxxx' to xxxx


class ThreadWorkerForWidgetEdit(threading.Thread):
  index = [0] #need unambigous name...

  def __init__(self, pipe, edit, color="Black", emit=True):
    super(ThreadWorkerForWidgetEdit, self).__init__()
    self._objectName = "ThreadWidgetEdit" + str(self.index)
    self.index[0] += 1
    #if verbose: print "ThreadWorkerForWidgetEdit.__init__",edit
    self.pipe = pipe
    self.edit = edit
    self.editObjectName = self.edit.objectName()
    self.color = color
    #http://stackoverflow.com/questions/323972/is-there-any-way-to-kill-a-thread-in-python
    self._stop = threading.Event()
    if verbose: print("DEBUG: init %s for %s"%(self._objectName,  self.editObjectName))
    #self.setDaemon(True) #may be serve: nokill on exit?
    self.emit = emit

  def objectName(self):
    return self._objectName
    
  def run(self):
    if verbose: print("DEBUG: run %s for %s"%(self._objectName,  self.editObjectName))
    self.worker()

  def worker(self):
    """
    http://stackoverflow.com/questions/2104779/qobject-qplaintextedit-multithreading-issues
    Although QObject is reentrant, the GUI classes, 
    notably QWidget and all its subclasses, are not reentrant.
    They can only be used from the main thread.
    i.e. use signals...
    warning 2017 : I never seen finishPopen on stderr don't know why...
    """
    verb = verbose
    windows_line_ending = '\r\n'
    linux_line_ending = '\n'
    while True:
      line = self.pipe.readline()
      try:
        line = line.decode("utf-8") # python3 stdout is b'...'
      except:
        print("DEBUG: decode problem ThreadWorkerForWidgetEdit '%s'" % line[:-1]) # minus "\n"
        line = str(line)

      line = line.replace(windows_line_ending, linux_line_ending)
      if len(line) < 10: # debug supposed line Popen finished
        if verb: print("DEBUG: Worker %s get '%s' as '%s'" % (self._objectName, line[:-1], toHexAll(line)))
      else:
        if verb: print("DEBUG: Worker %s get '%s'" % (self._objectName, line[:-1]))

      if line == "":
        if verb: print("DEBUG: Worker Popen finished %s for %s" % (self._objectName,  self.editObjectName))
        #sometimes error pyqtSignal must be bound to a QObject, not 'QTextEditForLog'
        #may be reentrance on self.edit delete?
        self.pipe.close() #add by ab alex berlin
        try:
          if self.emit: self.edit.finishPopen.emit(self._objectName)
          break
        except:
          if verb: print("DEBUG: Worker problem self.edit.finishPopen.emit", type(self.edit))
          break
      if self.stopped(): 
        if verb: print("DEBUG: Worker stopped %s for %s" % (self._objectName,  self.editObjectName))
        break
      else:
        try:
          self.edit.trigger.emit(line, self.color)
        except:
          print("WARNING: Worker problem self.edit.trigger.emit", type(self.edit))
          break

  def stop(self):
    if verbose: print("WARNING: ThreadWorkerForWidgetEdit have to stop %s for %s"%(self._objectName,  self.editObjectName))
    self._stop.set()

  def stopped(self):
    return self._stop.isSet()
