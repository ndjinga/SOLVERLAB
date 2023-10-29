#!/usr/bin/env python
# -*- coding: utf-8 -*-

# %% LICENSE_SALOME_CEA_BEGIN
# see PACKAGESPY/LICENCE file
# %% LICENSE_END


import os
import fnmatch
from time import sleep
from PyQt5 import QtGui, QtWidgets
from PyQt5 import QtCore

from salomepy.onceQApplication import OnceQApplication
import xyzpy.loggingXyz as LOG
from salomepy.qTextEditForLog import QTextEditForLog
import xyzpy.guiXyz.FileSystemXyz as FISY


logger = LOG.getLogger()
verbose = False
verboseEvent = verbose

"""
widget for a centralWidget Qmainwindow: multiple tab with editors for files and logs
"""

########################################################################################
class QTabMultipleTextCentral(QtWidgets.QTabWidget):
  
  #-1 as useless as not created
  #user use inheritage for modify TAB_... order or appearence tabs
  TAB_LOG_CMD = 0
  TAB_XML = 1
  TAB_OTHERTEXTEDIT = 2
  TAB_FILESYSTEM = -1 
  TAB_LOG_CAST3M = -1 #3 #-1 as useless as not created
  TAB_DGIBI = -1 #4
  TAB_PY = -1 #5
  TAB_IMAGEEDIT = -1 #6
  TAB_CSVEDIT = -1 #7

  index = [0]
  
  def __init__(self):
    super(QTabMultipleTextCentral, self).__init__()
    self.setObjectName("QTabMultipleTextCentral"+str(self.index))
    self.setWindowTitle(self.objectName())
    self.index[0] += 1
    
    self.activeTabs = {}
    self._init_tabs_standart()
    self._init_tabs_inherited()

    for k in sorted(self.activeTabs.keys()):
      wid = self.activeTabs[k]
      self.insertTab(k, wid, wid.tabName)
      setattr(self, wid.attName, wid)

    self.logSalomeWidget = self.logCMDWidget #salome is CMD appli
    self.resize(800,500)
    if verbose: 
      res = [(k, w.tabName) for k, w in list(self.activeTabs.items())]
      print("QTabMultipleTextCentral tabs: %s" % res)

  def _init_tabs_inherited(self):
    """virtual method for inherited additional tabs""" 
    pass
   
  def _init_tabs_standart(self):
   
    if self.TAB_LOG_CAST3M != -1:
      wid = QTextEditForLog()
      wid.saveFileExt = ".log"
      wid.tabName = "Log Castem" #as foreign clear short label name
      wid.attName = "logCastemWidget" #as immmutable name
      self.activeTabs[self.TAB_LOG_CAST3M] = wid
    
    if self.TAB_DGIBI != -1: 
      wid = QTextEditForLog()
      wid.saveFileExt = ".dgibi"
      wid.tabName = "File Dgibi"
      wid.attName = "dgibiWidget"
      self.activeTabs[self.TAB_DGIBI] = wid

    if self.TAB_PY != -1: 
      wid = QTextEditForLog()
      wid.saveFileExt = ".py"
      wid.tabName = "File py"
      wid.attName = "pyWidget"
      self.activeTabs[self.TAB_PY] = wid

    if self.TAB_XML != -1: 
      wid = QTextEditForLog()
      wid.saveFileExt = ".xml"
      wid.tabName = "File xml"
      wid.attName = "xmlWidget"
      self.activeTabs[self.TAB_XML] = wid

    if self.TAB_LOG_CMD != -1: 
      wid = QTextEditForLog()
      wid.saveFileExt = ".log"
      wid.tabName = "Log Run Code"
      wid.attName = "logCMDWidget"
      self.activeTabs[self.TAB_LOG_CMD] = wid

    if self.TAB_OTHERTEXTEDIT != -1: 
      wid = QTextEditForLog()
      wid.saveFileExt = ".tmp"
      wid.tabName = "Other File"
      wid.attName = "otherTexteditWidget"
      self.activeTabs[self.TAB_OTHERTEXTEDIT] = wid

    if self.TAB_FILESYSTEM != -1: 
      wid = FISY.FileSystemModelViewerWidget()
      #wid.saveFileExt = ".tmp"
      wid.tabName = "Explore Dir"
      wid.attName = "exploreDirWidget"
      self.activeTabs[self.TAB_FILESYSTEM] = wid
    
  def getTabs(self):
    return [self.activeTabs[k] for k in sorted(self.activeTabs.keys())]

  def getTabByName(self, aName):
    """
    name as attribute name or name as tab label
    wid.tabName = "Explore Dir"      #tab label cf tr()
    wid.attName = "exploreDirWidget" #immutable attribute name
    """
    tabs = self.getTabs()
    for wid in tabs:
       if aName == wid.tabName: return wid
       if aName == wid.attName: return wid
    return None

  def closeEvent(self, event):
    if verboseEvent: print(self.objectName()+".closeEvent")
    for k in sorted(self.activeTabs.keys()):
      wid = self.activeTabs[k]
      wid.close()
    return super(QTabMultipleTextCentral, self).closeEvent(event)
    
  def showTab(self, index):
    if type(index) == int:
      if index in list(self.activeTabs.keys()):
        self.setCurrentIndex(index)
      else:
        logger.warning("index %i unknown" % index)
    if type(index) == str:
      if index in list(self.activeTabs.keys()):
        self.setCurrentIndex(index)
      else:
        logger.warning("index %i unknown" % index)

  def showTabByName(self, aName):
    indexFound = None
    for index, wid in list(self.activeTabs.items()):
       if aName == wid.tabName or aName == wid.attName: 
         indexFound = index
         break
    if indexFound != None:
      self.setCurrentIndex(indexFound)
    else:
      logger.warning("tab '%s' unknown" % aName)

  def showLogCastemWidget(self):
    self.showTab(self.TAB_LOG_CAST3M)

  def showDgibiWidget(self):
    self.showTab(self.TAB_DGIBI)

  def showPyWidget(self):
    self.showTab(self.TAB_PY)

  def showXmlWidget(self):
    self.showTab(self.TAB_XML)

  def showLogSalomeWidget(self):
    self.showLogCMDWidget()

  def showLogCMDWidget(self):
    self.showTab(self.TAB_LOG_CMD)

  def showOtherTexteditWidget(self):
    self.showTab(self.TAB_OTHERTEXTEDIT)

  def showImageEditWidget(self):
    self.showTab(self.TAB_IMAGEEDIT)

  def showCsvEditWidget(self):
    self.showTab(self.TAB_CSVEDIT)

  def showExploreDirWidget(self):
    self.showTab(self.TAB_FILESYSTEM)

  """def hideLogSalomeWidget(self):
    self.hideLogCMDWidget()

  def hideOtherTexteditWidget(self):
    self.removeTab(self.TAB_OTHERTEXTEDIT)

  def hideXmlWidget(self):
    self.removeTab(self.TAB_XML)"""

  def quickEditFiles(self, files, mdump_all=None):
    logger.debug("Files %s",files)
    for i in files:
      """
      if fnmatch.fnmatch(i, "*.dgibi"):
        self.dgibiWidget.openFile(i)
        self.showDgibiWidget()
        continue
      if fnmatch.fnmatch(i, "*.py"):
        self.patternWidget.openFile(i)
        self.showPatternWidget()
        continue
      """
      if fnmatch.fnmatch(i, "*.xml"):
        self.xmlWidget.openFile(i)
        self.showXmlWidget()
        continue
      if fnmatch.fnmatch(i, "*.med"):
        # obsolete
        # f = open('/tmp/echo111.tmp','w')
        # f.write('1\n1\n1\n0\n')
        # f.close()
        user = os.getenv("USER", "unknownUser")
        basename = "%s_%s_Mdumped" % (os.path.basename(i), user)
        tmpFile = "/tmp/" + basename
        if mdump_all is None:
          MDUMP_ALL = os.getenv("MDUMP_ALL", None)
        else:
          MDUMP_ALL = ""
        if str(MDUMP_ALL) in "y Y yes YES 1 True".split():
          structure = ""
        else:
          structure = "--structure"  # sans afficher les données volumineuses
        """
        mdump [--structure] monfichier.med [ NODALE|DESCENDANTE NO_INTERLACE|FULL_INTERLACE|LECTURE_EN_TETE_SEULEMENT N°MAILLAGE|0 pour tous ]
        	--structure               : Lis l'ensemble des données sans afficher les données volumineuses
            NODALE                    : Scrute la connectivité nodale      des maillages
            DESCENDANTE               : Scrute la connectivité descendante des maillages
            FULL_INTERLACE            : Affiche les connectivités en mode      entrelacé x1y1x2y2
            NO_INTERLACE              : Affiche les connectivités en mode  non entrelacé x1x2y1y2
            LECTURE_EN_TETE_SEULEMENT : Affiche uniquement les entêtes, désactive la lecture et l'affichage des données volumineuses
            N°MAILLAGE ==  i          : Affiche le maillage n°i et ses champs associés
            N°MAILLAGE ==  0          : Affiche l'ensemble des maillages et leurs champs associés
            N°MAILLAGE == -1          : Affiche l'ensemble des champs qu'ils soient associés ou non à un maillage local
        """
        mdump = os.path.expandvars("${MEDFILE_ROOT_DIR}/bin/mdump")
        cmd = "%s %s %s NODALE NO_INTERLACE 0 > %s" % (mdump, structure, i, tmpFile)
        logger.info("mdump command:\n%s" % cmd)
        self.otherTexteditWidget.launchIntoPopen(cmd)
        sleep(1) # seems to be util
        self.otherTexteditWidget.openFile(tmpFile)
        self.showOtherTexteditWidget()
        continue
      """
      if fnmatch.fnmatch(i, "*.png"):
        self.showImageEditWidget()
        self.imageEditWidget.openFile(i)
        #self.showOtherTexteditWidget()
        continue
      if fnmatch.fnmatch(i, "*.help.html"): #pattern help files
        self.showImageEditWidget()
        self.imageEditWidget.openFile(i)
        #self.showOtherTexteditWidget()
        continue
      """
      self.otherTexteditWidget.openFile(i)
      self.showOtherTexteditWidget()

if __name__ == '__main__':
  app = OnceQApplication([''])
  fen = QTabMultipleTextCentral()
  #fen.logSalomeWidget.launchIntoPopen("pwd ; ls -alt *py ; someErrorExampleInRed")
  fen.showLogSalomeWidget()
  fen.show()
  print("finish, here comes app.exec_()")
  app.exec_()
  print("===Bye World===")

