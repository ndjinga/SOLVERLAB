#!/usr/bin/env python
# -*- coding: utf-8 -*-

# %% LICENSE_SALOME_CEA_BEGIN
# see PACKAGESPY/LICENCE file
# %% LICENSE_END


"""
temporary same as matplotbibWindowToolbar
later may be tabwidget with multiple MatplotlibWindowToolbar(s)
"""

from matplotlibpy.matplotlibWindowToolbar import MatplotlibWindowToolbar

verboseEvent = True

class QMainWindowForPlt(MatplotlibWindowToolbar):

  def receiveRequestToView(self, strXmlRequest):
    if verboseEvent: 
      print("%s %s receiveRequest virtual" % (self.__class.__name__, self.objectName()))
    return True
