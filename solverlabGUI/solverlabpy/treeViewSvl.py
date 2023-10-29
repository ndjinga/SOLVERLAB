#!/usr/bin/env python
# -*- coding: utf-8 -*-

# %% LICENSE_SALOME_CEA_BEGIN
# Copyright (C) 2008-2023  CEA
#
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 2.1 of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with this library; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
#
# See http://www.salome-platform.org or email : webmaster.salome@opencascade.com
# %% LICENSE_END


import os
import sys
import pprint as PP
from PyQt5 import QtGui, QtCore, QtWidgets as QTW

from xyzpy.guiXyz.treeXmlXyz import TreeXmlXyz, TreeXmlXyzItem
import xyzpy.loggingXyz as LOG
import solverlabpy.abcdExpression as ABCD
import solverlabpy.configSvl as CFGSVL

logger = LOG.getLogger()
verbose = False
verboseEvent = False

"""
cosmetic stuff for treeView Solverlab
"""

FORMATS_TREEVIEW = {
  "FMT_FLOAT": "{0:13.5g}",  # or "{0:13.5e}"
  "FMT_INT": "{0:>13}",
}


####################################################
class TreeViewSvl(TreeXmlXyz):
  class COLS:
    labels = ['Name', 'Value', 'Attributes']
    Tag = 0
    Text = 1
    Attributes = 2

  def __init__(self, parent=None):
    super(TreeViewSvl, self).__init__(parent)

    self.setHeaderLabels(self.COLS.labels)
    """
    #set in app.setFont
    self.setFont(QtGui.QFont("Monospace", 9))
    font = QtGui.QFont(self.font()) #copy: self.font() is const
    font.setFamily("Monospace")
    font.setPointSize(9)
    self.setFont(font)
    """
    self.setAlternatingRowColors(True)
    pal = self.palette()
    config = CFGSVL.getMainConfigCatchAll()
    if verbose: print("current configuration:\n%s" % config)
    colb = CFGSVL.toListInt(config.MainWindow.color_treeview_base)
    colt = CFGSVL.toListInt(config.MainWindow.color_treeview_text)
    pal.setColor(pal.Base, QtGui.QColor(*colb))
    pal.setColor(pal.Text, QtGui.QColor(*colt))
    self.setPalette(pal)
    # TODO remove that later self._TreeXmlXyzItemClass = TreeXmlXyzItemSvl
    # default formats for treeViewSvl
    self.formats_treeview = FORMATS_TREEVIEW # e to g format



