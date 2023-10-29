#!/usr/bin/env python
# -*- coding: utf-8 -*-

# %% LICENSE_SALOME_CEA_BEGIN
# see PACKAGESPY/LICENCE file
# %% LICENSE_END

import os
import sys
import traceback

#from PyQt5.QtGui import *
#from PyQt5.QtCore import *
import PyQt5.QtWidgets as QTW

import xyzpy.loggingXyz as LOG

logger = LOG.getLogger()

verbose = False

###############################################################
class MessageBoxDialog(QTW.QMessageBox):

  def resizeEvent(self, Event):
    # logger.info("MessageBoxDialog.resizeEvent")
    QTW.QMessageBox.resizeEvent(self, Event)
    self.setFixedWidth(700)


###############################################################
def getMessageBoxDialog(parent=None, mess="", level="ERROR", details=None, exception=None, icon="Warning"):
  wid = MessageBoxDialog(parent)
  wid.setText(level)
  wid.setIcon(wid.Warning)
  allmess = mess
  if icon is not None:
    if type(icon) is str:
      iconl = icon.lower()
      if iconl == "warning": iconc = wid.Warning
    else:
      iconc = icon
    wid.setIcon(iconc)
  if exception is not None:
    allmess = mess + '\n' + str(exception)
  wid.setInformativeText(allmess + "\n") # space ...
  if details is not None:
    wid.setDetailedText(details)
  else:
    trace = traceback.format_exc() #better explicit verbose problem
    wid.setDetailedText(trace)
  return wid
