#!/usr/bin/env python
# -*- coding: utf-8 -*-

# %% LICENSE_SALOME_CEA_BEGIN
# see PACKAGESPY/LICENCE file
# %% LICENSE_END

"""
as QtWidgets.QFileDialog() with setSidebarUrls with environment variables

>>> from widgetpy.salomeQFileDialog import SalomeQFileDialog
>>> aDialog=SalomeQFileDialog()
>>> aFile=aDialog.browseFileDialog("aTitle", "aFileOrDirIni", "(*)", ['HOME','MATIXWORKDIR'])
>>> aDir=aDialog.browseDirDialog("aTitle", "aFileOrDirIni", ['MATIXWORKDIR'], showDirsOnly=False)
"""

import os
from PyQt5 import QtGui,QtCore, QtWidgets


class SalomeQFileDialog(QtWidgets.QFileDialog):
  
  def __init__(self, parent=None):
    super(SalomeQFileDialog, self).__init__(parent)
    self._envvars = ['HOME','MATIX_WORKDIR']
    self._title = 'Select file'
    self._filter = '(*)' #"(*.xml *.XML)")
    self._nameFile = ""
    self._nameDir = ""
    self._verbose = False
    
  def resetSidebarUrls(self, aEnvVars=['HOME','MATIX_WORKDIR']):
    """
    force SidebarUrls because sometimes too much items staying alive
    TODO problem set SidebarUrls have to do 2 times, don't know why see self._verbose
    """
    urlsIni = [i for i in self.sidebarUrls()]
    if self._verbose: print("urlsIni",urlsIni)
    urls = [QtCore.QUrl('file:')] #unremovable
    for var in aEnvVars:
      aDir = os.getenv(var)
      if aDir != None:
        oneUrl = QtCore.QUrl.fromLocalFile(aDir)
        if oneUrl not in urls: urls.append(oneUrl)
    verb = False #self._verbose
    if verb:
      print("--reset urlssidebar before",self.sidebarUrls())
      print("--reset urlssidebar new set",urls)
    self.setSidebarUrls(urls)
    if verb: 
      print("--reset urlssidebar after",self.sidebarUrls())
    self.update()

  def browseFileDialog(self, aTitle=None, aFileOrDir=None, aFilter=None, aEnvVars=None):
    import salomepy.xsalomesession as XSS

    if self._verbose: print("browseFileDialog parent %s" % self.parent())

    aTitlec = self._title
    if aTitle != None: aTitlec=str(aTitle)
    
    aFileOrDirc = ""
    if aFileOrDir != None: aFileOrDirc=str(aFileOrDir)
    
    aFilterc = self._filter
    if aFilter != None: aFilterc=str(aFilter)
    
    aEnvVarsc = [var for var in self._envvars] #copy
    if aEnvVars != None:
      for var in aEnvVars:
        aEnvVarsc.append(var)

    try:
      desktop = XSS.getDesktop()
    except:
      desktop = self.parent() # None  #if not context salome
    urlsIni = [QtCore.QUrl(i) for i in self.sidebarUrls()] #copy
    if self._verbose: print("urlsIni",urlsIni)
    urls = [QtCore.QUrl(i) for i in urlsIni] #copy

    for var in aEnvVarsc:
      aDir = os.getenv(var)
      if aDir != None:
        oneUrl = QtCore.QUrl.fromLocalFile(aDir)
        if oneUrl not in urls: urls.append(oneUrl)

    if aFileOrDirc != "":
      if os.path.isdir(aFileOrDirc):
        aDir = aFileOrDirc
      else:
        aDir, _ = os.path.split(aFileOrDirc)
      oneUrl = QtCore.QUrl.fromLocalFile(aDir)
      if oneUrl not in urls: urls.append(oneUrl)
    
    if self._verbose: print("browseFileDialog urls %s" % urls)

    self.setSidebarUrls(urls)
    if self._verbose: print("urlssidebar before %s" % self.sidebarUrls())
    
    #https://stackoverflow.com/questions/43509220/qtwidgets-qfiledialog-getopenfilename-returns-a-tuple
    res = self.getOpenFileName(desktop, aTitlec, str(aDir), aFilterc)
    nameFile = str(res[0])

    if self._verbose: print("open file %s" % str(res))

    self._nameFile = nameFile
    if nameFile != "":
      realPath = os.path.realpath(nameFile)
      self._nameFile = realPath
      if not os.path.isfile(realPath):
        QtWidgets.QMessageBox.warning(desktop, "warning", "%s: not a file \n'%s'" % (aTitlec, realPath))
        self._nameFile = ""
    #self.setSidebarUrls(urlsIni)
    if self._verbose: print("urlssidebar after",self.sidebarUrls())
    return self._nameFile

  def browseDirDialog(self, aTitle=None, aFileOrDir=None, aEnvVars=None, showDirsOnly=False):
    import salomepy.xsalomesession as XSS
    
    aTitlec = self._title
    if aTitle != None: aTitlec=str(aTitle)
    
    aFileOrDirc = ""
    if aFileOrDir != None: aFileOrDirc=str(aFileOrDir)
    
    #aFilterc = self._filter #no needs
    #if aFilter != None: aFilterc=str(aFilter)
    
    aEnvVarsc = [var for var in self._envvars] #copy
    if aEnvVars != None:
      for var in aEnvVars:
        aEnvVarsc.append(var)
    
    showDirsOnlyc=self.ShowDirsOnly
    if not showDirsOnly:
      showDirsOnlyc=self.Options(0) #http://pyqt.sourceforge.net/Docs/PyQt4/qfiledialog.html#Option-enum

    desktop = XSS.getDesktop()
    urlsIni = [QtCore.QUrl(i) for i in self.sidebarUrls()] #copy
    urls = [QtCore.QUrl(i) for i in urlsIni] #copy

    for var in aEnvVarsc:
      aDir = os.getenv(var)
      if aDir != None:
        oneUrl = QtCore.QUrl.fromLocalFile(aDir)
        if oneUrl not in urls: urls.append(oneUrl)

    if aFileOrDirc != "":
      if os.path.isdir(aFileOrDirc):
        aDir = aFileOrDirc
      else:
        aDir, _ = os.path.split(aFileOrDirc)
      oneUrl = QtCore.QUrl.fromLocalFile(aDir)
      if oneUrl not in urls: urls.append(oneUrl)
    
    if self._verbose: print("urls",urls)
    self.setSidebarUrls(urls)
    if self._verbose: print("urlssidebar before",self.sidebarUrls())
    
    nameDir = str(self.getExistingDirectory(desktop, aTitlec, str(aDir), showDirsOnlyc))
    self._nameDir = nameDir
    if nameDir != "":
      realPath = os.path.realpath(nameDir)
      self._nameDir = realPath
      if not os.path.isdir(realPath):
        QtWidgets.QMessageBox.warning(desktop, "warning", "%s: not a directory \n'%s'" % (aTitlec, realPath))
        self._nameDir = ""
    #self.setSidebarUrls(urlsIni)
    if self._verbose: print("urlssidebar after",self.sidebarUrls())
    return self._nameDir


