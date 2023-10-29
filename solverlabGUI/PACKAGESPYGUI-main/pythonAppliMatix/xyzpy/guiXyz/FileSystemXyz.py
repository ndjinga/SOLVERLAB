#!/usr/bin/env python
# -*- coding: utf-8 -*-

# %% LICENSE_SALOME_CEA_BEGIN
# see PACKAGESPY/LICENCE file
# %% LICENSE_END


"""
widget directories and files explorer and selecting functionalities with 
- syntax highlight viewer ascii
- viewer/browser html
"""

import os
import sys
import fnmatch
import platform

from PyQt5 import QtCore, QtGui, QtWidgets
  
import subprocess as SP
import mimetypes
import pprint as PP

from salomepy.onceQApplication import OnceQApplication
import salomepy.iconsUser as IUSR
import xyzpy.loggingXyz as LOG
import xyzpy.utilsXyz as UXYZ #common procedures
import widgetpy.messageBoxDialog as MBD


logger = LOG.getLogger()

# http://doc.qt.io/qt-5/qtwebenginewidgets-qtwebkitportingguide.html
# The WebKit is indeed not supported anymore since Qt5.6
# have to set PyQt5.11 ?
# PyQt5.QtCore.QT_VERSION_STR -> 5.6.2 for now
# PyQt5.QtWebKit.qWebKitVersion() -> 538.1
try:
  # avoid message under windows 'Qt WebEngine seems ... from a plugin'
  # from PyQt5.QtWebEngineWidgets import QWebEngineView
  logger.info("PyQt5.QtWebEngineWidgets.QWebEngineView is not used (may be available)")
  QWebEngineView = None # not used
  # from PyQt5 import QtWebKit
  # from PyQt5.QtWebKit import QWebView
except:
  QWebEngineView = None
  logger.critical("No PyQt5.QtWebEngineWidgets.QWebEngineView available")


verbose = False
verboseView = False  # view user event
verboseEvent = False  # very verbose

testDir = os.path.realpath(
          os.path.split(
          os.path.realpath(__file__))[0])

_themesDir = "/usr/share/highlight/themes"

# risky on big files
mimetypes.add_type('application/uranie', '.dat') # unknown
mimetypes.add_type('notasci/corteo', '.asp') # binary data from corteo for iradina
mimetypes.add_type('notasci/bz2', '.bz2') # binary data from zip
mimetypes.add_type('notasci/7z', '.7z') # binary data from zip
mimetypes.add_type('text/iradina', '.cfg') # asci config file for iradina
mimetypes.add_type('text/iradina', '.in') # asci config file for iradina
mimetypes.add_type('text/iradina', '.ini') # asci config file for iradina
mimetypes.add_type('text/iradina', '.log') # asci config file for iradina

############################################
def isTextFileMimetypes(path):
  """if file acceptable for an ascii text editor"""
  base = os.path.basename(path)
  mime = mimetypes.guess_type(base)
  logger.info("mimetype: %s %s" % (mime, path))

  include_pattern = "ira.* *.dgibi *.trace *.bash *.log".split()
  for i in include_pattern:
    if fnmatch.fnmatch(base, i):
      logger.debug("file %s as text from pattern '%s'" % (base, i))
      return True

  # ??? linux mimetypes.add_type('notasci/bz2', '.bz2') do not work
  exclude_ext = [".bz2"]
  ext = os.path.splitext(path)[1]
  if ext in exclude_ext:
    return False

  # '.so' as 'application/octet-stream'  '.sh' as 'application/x-sh'
  exclude = ["image", "x-python-code", "octet-stream", "notasci", "zip"]
  mim0 = mime[0]
  if mim0 == None: 
    logger.warning("unknown mimetype file %s, try (risky) as text." % base)
    return True #unknown mime as text editable (risky)
  for ex in exclude:
    if ex in mime[0]: return False
  return True

############################################
def isImageFileMimetypes(path):
  """is file image as tiff png jpeg etc..."""
  mime = mimetypes.guess_type(path)
  if verbose: logger.info("mimetype: %s %s " % (mime, os.path.basename(path)))
  mim0 = mime[0]
  if mim0 == None: return False #unknown mime as not image
  if "image" in mime[0]: return True
  return False

############################################
def isHtmlFileMimetypes(path):
  """is file html for direct display in web browser"""
  mime = mimetypes.guess_type(path)
  if verbose: logger.info("mimetype: %s %s " % (mime, os.path.basename(path)))
  mim0 = mime[0]
  if mim0 == None: return False #unknown mime as not html
  if "html" in mime[0]: return True
  return False

############################################
def isWebFileMimetypes(path):
  """is file acceptable for a web browser"""
  if isHtmlFileMimetypes(path): return True
  if isImageFileMimetypes(path): return True
  return False

############################################
def isLibraryFile(path):
  if os.path.splitext(path)[1] in [".so"]: return True
  return False

############################################
def isIradinaOutputFile(nameFile):
  """
  assert iradina output files from base name of file
  if prefix is 'test.' 'ira.' or else...
  as 'ira.etc.etc' or 'test.etc.etc'
  """
  prefix = os.path.basename(nameFile).split(".")[0]  # ira.etc.etc or test.etc.etc
  if prefix in "test ira".split():
    return True
  return False

############################################
def isIradinaInFile(nameFile):
  """
  assert iradina inputput files from base name of file
  if extension is '.cfg' '.in'
  as 'Configuration.in'
  """
  ext = os.path.splitext(nameFile)[1]  # get extension file name (as '.in')
  if ext in ".cfg .in .log".split():
    return True
  return False

############################################
def isTextPrefixFile(nameFile):
  """
  assert text output files from base name of file
  if prefix is 'test.' 'ira.' or else...
  as 'ira.etc.etc' or 'test.etc.etc'
  """
  prefix = os.path.basename(nameFile).split(".")[0]  # ira.etc.etc or test.etc.etc
  if prefix in "test ira out".split():
    return True
  return False

############################################
def isIniExtensionFile(nameFile):
  """
  assert text files from base name of file
  if extension is '.cfg' '.in'
  as 'Configuration.in'
  """
  ext = os.path.splitext(nameFile)[1]  # get extension file name (as '.in')
  if ext in ".cfg .in .dgibi .trace .log".split():
    return True
  return False

############################################
def isTextExtensionFile(nameFile):
  """
  assert text files
  """
  ext = os.path.splitext(nameFile)[1]  # get extension file name (as '.in')
  if ext in ".cfg .in .log .dgibi .bash .trace".split():
    return True
  return False

############################################
def _createAction(anObject, Name, Shortcut, ToolTip, Call, Icon=None, Enable=True):
  """utility for classes"""
  if verbose: print("FileSystem._createAction for %s '%s" % (anObject.__class__.__name__, Name))
  action = QtWidgets.QAction(Name, anObject)
  action.setIconVisibleInMenu(Icon != None)
  if Shortcut != None: 
    action.setShortcut(anObject.prefixShortcut+Shortcut)
  action.setToolTip(ToolTip)
  if Icon != None:
    action.setIcon(IUSR.getIconFromName(Icon))
  #anObject.connect(action, QtCore.SIGNAL('triggered()'), Call) #old design
  action.setEnabled(Enable)
  action.triggered.connect(Call)
  return action


############################################
class DirView(QtWidgets.QTreeView):

  SetRootPathSignal = QtCore.pyqtSignal(str)

  def __init__(self, *args, **kwargs):
    super(DirView, self).__init__(*args, **kwargs)
    self._createActions()
    self._createContextMenus()

  def getWindowForLog(self):
    """example not serious for debug"""
    if True: #try:
      ok = False
      parent = self
      for i in range(20):
        # print "DirView.getWindowForLog", parent.__class__.__name__
        if parent == None: return None
        if "QMainWindowForLog" in parent.__class__.__name__: return parent
        parent = parent.parent()
    return None

  def event(self, event):
    if verboseEvent: logger.info("DirView.event %s" % event)
    return super(DirView, self).event(event)

  def _createActions(self):
    self.expandActions = [
      _createAction(self, 'Expand all', None, 'expand all items', self.expandAll, 'expand'),
      # _createAction(self, 'display ROOT(CERN) directory', None, 'display ROOT(CERN) directory contents', self.emitDisplayRoot, 'root'),
      _createAction(self, 'display parent root directory', None, 'set root directory', self.emitDisplayParentRoot, 'insertitemup'),
      _createAction(self, 'browse new root directory', None, 'set new root directory', self.browseRootDirectory, 'browsefile'),
      # TODO _createAction(self, 'search in files', None, 'search string in files (recursive)', self.findInFiles, 'search'),
    ]

  def _createContextMenus(self):
    """each column have a different menu, virtual menus as example for inheritage"""
    menus = {}
    
    #outside column event
    menu = QtWidgets.QMenu("General", self)
    for action in self.expandActions: 
      menu.addAction(action)
    menus[-1] = menu
    
    #column 0 event
    #menu = QtWidgets.QMenu("Menu_0", self)
    #for action in #TODO: 
    #  menu.addAction(action)
    menus[0] = menu
    
    #column 1 event
    #menu = QtWidgets.QMenu("Menu_1", self)
    menus[1] = menu
    self.contextMenus = menus
    return

  def contextMenuEvent(self, event):
    """
    each column have a different menu, reals menus from user inheritage
    from user implementation of _createContextMenus
    """
    index = self.selectionModel().currentIndex()
    row = index.row()
    column = index.column()
    if verboseView: logger.info("DirView.contextMenuEvent row %i column %i" % (row,column))
    
    pos = QtGui.QCursor.pos()
    theMenu = self.contextMenus[column]
    #theMenu.exec_(self.mapToGlobal(QtCore.QPoint(10,10)))
    if verboseView: logger.info("DirView parent widget '%s'" % self.mainWidget._name)
    theMenu.exec_(pos)
    return

  def emitDisplayRoot(self):
    self.clicked.emit(None) #as root

  def emitDisplayParentRoot(self):
    self.SetRootPathSignal.emit("UP") #as Parent of root

  def browseRootDirectory(self):
    dial = QtWidgets.QFileDialog()
    directory = dial.getExistingDirectory(self,
                "choose new root directory", "/", dial.ShowDirsOnly)
    if directory:
      self.SetRootPathSignal.emit(str(directory))

  def refreshFileView(self):
    selection = self.selectionModel()
    if selection.hasSelection():
      index = selection.currentIndex()
      self.clicked.emit(index) #will refresh fileView
    else:
      # https://stackoverflow.com/questions/12050397/pyqt-signal-emit-with-object-instance-or-none
      self.clicked.emit(type(None)) #as root

  def findInFiles(self):
    """
    quick linux recursive find string in files from root
    output console, for the moment
    """
    dial = QtWidgets.QInputDialog
    currentDir = self.model().rootPath() 
    # print "search in files in ",currentDir
    messageIni = "type here your mandatory search string"
    message, ok = dial.getText(self, "search recursive", 
                    "in %s" % currentDir, 
                    QtWidgets.QLineEdit.Normal, messageIni)
    if not ok: return False
    if message == messageIni: return False
    if message == "": #find all files
      cmd = "find . -type f"
    else:
      cmd = "cd %s ; fgrep -RHn -m 1 '%s'" % (currentDir, message)
    # rc = UXYZ.Popen(cmd, logger=logger)
    self.getWindowForLog().launchCmdIntoPopen(cmd)
    return 

############################################
class FileView(QtWidgets.QTreeView):
  
  helpSignal = QtCore.pyqtSignal()
  
  def __init__(self, *args, **kwargs):
    super(FileView, self).__init__(*args, **kwargs)
    #self.helpSignal.connect(self.on_help)
    #self.setToolTip("FileView")
    self._createActions()
    self._createContextMenus()
 
  def event(self, event):
    # if verboseEvent: print("FileView.event",event)
    if type(event) == QtGui.QHelpEvent: 
      # print "FileView help event", event.type()
      self.helpSignal.emit()
    return super(FileView, self).event(event)

  def _createActions(self):
    editor = UXYZ.getEditor(basename=True)
    browser = UXYZ.getBrowser(basename=True)
    if verboseView: logger.info("FileView._createActions")
    self.fileActions = [
      _createAction(self, '%s edit' % editor, None, 'edit by %s' % editor, self.editBy, 'edit'),
      _createAction(self, '%s view' % browser, None, 'view by %s' % browser, self.viewBy, 'webbrowser'),
      _createAction(self, 'display all files (*)', None, 'display all files', self.viewAll, None),
    ]
    self.plotActions = [
      _createAction(self, 'matplotlib plot', None, 'plot contents of file', self.emitPlotFile, 'plot'),
    ]
    for i in self.plotActions:
      self.fileActions.append(i)

  def _createContextMenus(self):
    """each column have a different menu, virtual menus as example for inheritage"""
    if verboseView: logger.info("FileView._createContextMenus")
    menus = {}
    
    # outside column event
    menu = QtWidgets.QMenu("General", self)
    for action in self.fileActions: 
      menu.addAction(action)
    menus[-1] = menu
    
    # column 0 event
    #menu = QtWidgets.QMenu("Menu_0", self)
    #for action in #TODO: 
    #  menu.addAction(action)
    menus[0] = menu
    
    # column 1 event
    #menu = QtWidgets.QMenu("Menu_1", self)
    menus[1] = menu
    self.contextMenus = menus
    return

  def contextMenuEvent(self, event):
    """
    each column have a different menu, reals menus from user inheritage
    from user implementation of _createContextMenus
    """
    if verboseView: logger.info("FileView.contextMenuEvent")
    index = self.selectionModel().currentIndex()
    row = index.row()
    column = index.column()
    tmp = self.model().fileInfo(index)
    # print("self.model().fileInfo(index) %s as %s" % (type(tmp), tmp))

    self.currentPath = str(self.model().fileInfo(index).absoluteFilePath())
    if verboseView:
      logger.info("FileView.contextMenuEvent row %i column %i\n  '%s'" % (row,column,self.currentPath))
    
    pos = QtGui.QCursor.pos()
    try:
      theMenu = self.contextMenus[column]
    except:
      # case unexpected column (> 1, see _createContextMenus)
      theMenu = self.contextMenus[1]

    # plot for iradina
    toEnable = os.path.basename(self.currentPath)[0:4] in "ira. test".split()
    for a in self.plotActions:
      a.setEnabled(toEnable)

    if verboseView: logger.info("FileView parent widget '%s'" % self.mainWidget._name)
    #theMenu.exec_(self.mapToGlobal(QtCore.QPoint(10,10)))
    theMenu.exec_(pos)
    return

  def emitPlotFile(self):
    if verboseView: logger.info("FileView.emitPlotFile %s" % self.currentPath)
    self.mainWidget.plotFileSignal.emit({"typePlot": "IradinaPlot", "filePlot": self.currentPath})

  def editBy(self):
    if verboseView: logger.info("FileView.editBy")
    if os.path.isdir(self.currentPath):
      QtWidgets.QMessageBox.warning(self,"warning","Problem is a directory:\n%s" % self.currentPath)
      return
    else:
      # may be if no echo then no return on communicate() under windows
      # windows7 'START /B' obsolescent?
      # import subprocess as SP
      # cmd = 'START /B %s %s' % (UXYZ.getEditor(), self.currentPath)
      cmd = '%s %s &' % (UXYZ.getEditor(), self.currentPath)
      # print("begin %s '%s'" % (__file__, cmd))
      # rc = UXYZ.Popen(cmd, logger=logger) # wait end
      res = SP.Popen(cmd, shell=True)
      # print("end '%s'" % str(res))

  def viewBy(self):
    if verboseView: logger.info("FileView.viewBy")
    cmd = "%s %s &" % (UXYZ.getBrowser(), self.currentPath)
    # rc = UXYZ.Popen(cmd, logger=logger) # wait end
    res = SP.Popen(cmd, shell=True)

  def viewAll(self):
    if verboseView: logger.info("FileView.viewAll")
    self.model().setNameFilters(["*"])
    if os.path.isfile(self.currentPath):
      aDir = os.path.dirname(os.path.realpath(self.currentPath))
    else:
      aDir = os.path.realpath(self.currentPath)
    cmd = "%s %s" % (UXYZ.getFileBrowser(), aDir)
    logger.info("viewAll file browser: '%s'" % cmd)
    process = SP.Popen(cmd, shell=True)

  """def on_help(self):
    index = self.selectionModel().currentIndex()
    row = self.selectionModel().currentIndex().row()
    column = self.selectionModel().currentIndex().column()
    model = self.model()
    #baseNode = self.selectedItems().getSelected[0]
    #for i in dir(self.selectionModel()): print "selectionModel",i
    #item = model.item(row,column)
    if verboseEvent: print "FileView.on_help index",row,column,model"""


############################################
class FileSystemModel(QtWidgets.QFileSystemModel):

  def __init__(self, *args, **kwargs):
    super(FileSystemModel, self).__init__(*args, **kwargs)
    #in case of use subset of selected files as FileItemDelegate (see delegate)
    self.selectedModel = None #case of use set of selected files

  def data(self, index, role):
    """to set tooltip"""
    if role == QtCore.Qt.ToolTipRole:
      # print "FileSystemModel.data", index, role
      model = index.model() #PyQt5.QtWidgets.QFileSystemModel
      path = model.fileInfo(index).absoluteFilePath()
      return path #yet QtCore.QVariant(path)
    return super(FileSystemModel, self).data(index, role)

  """do not work do not know why
  def setData(self, index, value, role):
    if True: #role == QtCore.Qt.ToolTipRole:
      print "FileSystemModel.setData", index, role
    return super(FileSystemModel, self).setData(index, value, role)"""


############################################
class FileSystemModelWidget(QtWidgets.QWidget):
  """
  a widget dialog with a view for dirs and another one for current files
  as http://doc.qt.io/qt-5/model-view-programming.htmlh
  """
  def __init__(self, *args):
    super(FileSystemModelWidget, self).__init__(*args)
    self.setWindowTitle("Files browser")
    self._name = "Files browser"
    layout = QtWidgets.QHBoxLayout()
    self.splitter = QtWidgets.QSplitter(QtCore.Qt.Horizontal)
    self.widget1 = QtWidgets.QWidget()
    self.widget2 = QtWidgets.QWidget()
    self.layout1 = QtWidgets.QVBoxLayout()
    self.layout2 = QtWidgets.QVBoxLayout()
    self.dirView = DirView()
    self.dirView.mainWidget = self
    self.fileView = FileView()
    self.fileView.mainWidget = self
    self.layout1.addWidget(self.dirView, 10) # 10 for dirView major occupation
    self.layout2.addWidget(self.fileView)
    self.widget1.setLayout(self.layout1)
    self.widget2.setLayout(self.layout2)
    self.splitter.addWidget(self.widget1)
    self.splitter.addWidget(self.widget2)
    layout.addWidget(self.splitter)
    self.setLayout(layout)
    self.selectedModel = None # a prori no subset of selected file

    self.theHeader=self.dirView.header()
    self.theHeader.setSectionsClickable(True)

    ##signals##
    self.theHeader.sectionClicked.connect(self.headerClikedEvent)
    self.theHeader.sectionDoubleClicked.connect(self.headerDoubleClickedEvent)

    # print "fileView.click",[i for i in dir(self.fileView) if "lick" in i]
    self.dirView.clicked.connect(self.on_dirView_clicked) #to refresh self.fileView
    self.dirView.SetRootPathSignal.connect(self.on_dirView_setRootPath) #to refresh self.fileView
    self.fileView.clicked.connect(self.on_fileView_clicked) #to refresh self.webView/editView
    self.fileView.doubleClicked.connect(self.on_fileView_doubleClicked) #for example for futur
    
  def setDirRootPath(self, rootPath, filters=[]):
    """
    rootPath is environ variables compliant as "$PACKAGESPY_ROOT_DIR/toto"
    filters as ["*.py"]
    """
    realRootPath = self.getRealPath(rootPath)
    if not os.path.isdir(realRootPath):
      raise Exception("Inexisting directory: '%s'" % realRootPath)
    self.originalRootPath = rootPath
    self.realRootPath = realRootPath 
    
    dirModel = FileSystemModel() #QtWidgets.QFileSystemModel()
    fileModel = FileSystemModel() #QtWidgets.QFileSystemModel()
    #a posteriori may be exist subset of selected
    fileModel.selectedModel= self.selectedModel

    dirModel.setReadOnly(True)
    fileModel.setReadOnly(True)
    
    rootDir = dirModel.setRootPath(self.realRootPath)
    rootFile = fileModel.setRootPath(self.realRootPath)
    
    dirModel.setFilter(QtCore.QDir.NoDotAndDotDot | QtCore.QDir.AllDirs)
    fileModel.setFilter(QtCore.QDir.NoDotAndDotDot | QtCore.QDir.Files)
    fileModel.setNameFilters(filters)
    fileModel.setNameFilterDisables(False) #else only disable (grey)
    
    self.dirView.setModel(dirModel)
    self.fileView.setModel(fileModel)

    self.dirView.setRootIndex(rootDir)
    self.fileView.setRootIndex(rootFile)

    self.dirModel = dirModel
    self.fileModel = fileModel

    self.dirView.hideColumn(1) #size
    self.dirView.hideColumn(2) #type
    self.dirView.hideColumn(3) #date modified
    
    self.fileView.header().setMinimumSectionSize(1)
    self.fileView.setColumnWidth(0, 200) #name
    self.fileView.setColumnWidth(1, 100) #size
    self.fileView.setColumnWidth(2, 1) #type
    self.fileView.setColumnWidth(3, 100) #date modified
    # self.fileView.hideColumn(2) #type
    
    self.fileView.setSortingEnabled(True)
    try:
      self.fileView.header().setSortIndicator(0, QtCore.Qt.AscendingOrder)
    except:
      #no header if self.fileView as QTableView
      pass

    #self.dirView.setRootIsDecorated(False) #useless
    #self.dirView.expandAll() #not working

    #self.setWindowTitle("Files browser (%s)" % self.realRootPath)
    self.setWindowTitle(self.realRootPath)

  def on_dirView_setRootPath(self, rootPath="UP"):
    if rootPath == "UP":
      path = os.path.realpath(os.path.split(self.realRootPath)[0])
    else:
      # print "on_dirView_setRootPath '%s' '%s'", type(rootPath), rootPath
      path = os.path.realpath(str(rootPath))
    # print "set parent RootPath", path
    self.setDirRootPath(path)

  def on_dirView_clicked(self, index):
    if index == None: #root path
      path = self.realRootPath
    else:
      path = self.dirModel.fileInfo(index).absoluteFilePath()
    if verboseEvent: logger.info("FileSystemModelWidget.on_dirView_clicked slot on\n  '%s'" % path)
    self.setRootFile(path)

  def on_fileView_clicked(self, index):
    path = self.fileModel.fileInfo(index).absoluteFilePath()
    if verboseEvent: logger.info("FileSystemModelWidget.on_fileView_clicked slot on\n'%s'" % path)
    if self.fileModel.isDir(index): #os.path.isdir(path): 
      self.setRootFile(path)
    else:
      #do nothing: edit/view if clicked see method inherited class FileSystemModelViewerWidget
      pass

  def on_fileView_doubleClicked(self, index):
    path = self.fileModel.fileInfo(index).absoluteFilePath()
    if verboseEvent: logger.info("FileSystemModelWidget.on_fileView_doubleClicked virtual slot on\n'%s'" % path)

  def event(self, event):
    # if verboseEvent: print("FileSystemModelWidget.event",event)
    return super(FileSystemModelWidget, self).event(event)

  def setRootFile(self, path):
    rootFile = self.fileModel.setRootPath(path)
    self.fileView.setRootIndex(rootFile)
    self.setWindowTitle(path)
    header = self.fileView.header()
    header.setStretchLastSection(False)
    # print("TODO setRootFile resizeColumnToContents")
    # self.fileView.resizeColumnToContents(0)
    # header.setSectionResizeMode(header.ResizeToContents)
    # header.setStretchLastSection(True)
    
  def headerClikedEvent(self, arg):
    if verboseEvent: logger.info("FileSystemModelWidget.headerClikedEvent at column %i" % arg)
    self.setRootFile(self.realRootPath)
    self.dirView.selectionModel().clearSelection()

  def headerDoubleClickedEvent(self, arg):
    if verboseEvent: logger.info("FileSystemModelWidget.headerDoubleClickedEvent at column %i" % arg)
    #self.setRootFile(self.realRootPath)
    #pos = QtGui.QCursor.pos()
    #theMenu = self.contextMenus[-1]
    #theMenu.exec_(self.mapToGlobal(QtCore.QPoint(10,10)))
    #theMenu.exec_(pos)

  def getRealPath(self, aPathWithEnvVar):
    """
    resolve env variable as ${HOME}/toto etc... 
    with os.path.expandvars interpretation of env var
    """
    res = os.path.expandvars(aPathWithEnvVar)  
    logger.debug("FileSystemModelWidget.getRealPath: '%s'->'%s'" % (aPathWithEnvVar, res))
    return res

  def toPathWithEnvVar(self, aPath):
    """
    returns a path with environ variable from original root path
    """
    aPathWithEnvVar = aPath.replace(self.realRootPath, self.originalRootPath)
    return aPathWithEnvVar

  
############################################
class FileItemDelegate(QtWidgets.QStyledItemDelegate):

  """to change paint in red if file is selected"""

  def paint(self, painter, option, index):

    modelItem = index.model() #PyQt5.QtWidgets.QFileSystemModel
    if modelItem.selectedModel == None:
      super(FileItemDelegate, self).paint(painter, option, index)
      return

    #modelItem.parent()
    #treeItem = treeWidget.itemFromIndex(index)
    #print "treeItem",treeItem,index.row(),index().column()
    
    path = modelItem.fileInfo(index).absoluteFilePath()
    
    toRed = modelItem.selectedModel.getStrList()
    textcolor = QtCore.Qt.black
    if path in toRed: textcolor = QtCore.Qt.red
   
    painter.save()

    """# set background color
    painter.setPen(QtGui.QPen(QtCore.Qt.NoPen))
    if option.state & QtWidgets.QStyle.State_Selected:
      painter.setBrush(QtGui.QBrush(QtCore.Qt.red))
    else:
      painter.setBrush(QtGui.QBrush(QtCore.Qt.white))
    painter.drawRect(option.rect)"""

    # set text color
    painter.setPen(QtGui.QPen(textcolor))
    value = index.data(QtCore.Qt.DisplayRole)
    painter.drawText(option.rect, QtCore.Qt.AlignLeft, value)
 
    painter.restore()
    return


############################################
class SelectedFilesListModel(QtCore.QStringListModel):
  #http://www.poketcode.com/en/pyqt4/item_views/index.html

  def data(self, index, role):
    """to set tooltip and display only basename of file"""
    #print "SelectedFilesListModel.data", index, role

    if role == QtCore.Qt.ToolTipRole:
      #yet QtCore.QVariant(path)
      return super(SelectedFilesListModel, self).data(index, QtCore.Qt.EditRole)

    if role == QtCore.Qt.DisplayRole:
      pathname = super(SelectedFilesListModel, self).data(index, QtCore.Qt.EditRole)
      #print "pathname '%s'" % pathname
      path, name = os.path.split(pathname)
      return QtCore.QVariant(name)

    return super(SelectedFilesListModel, self).data(index, role)

  def getStrList(self):
    """
    as list(str), copy NOT const
    http://doc.qt.io/qt-4.8/qstringlistmodel.html#stringList
    """
    return [str(i) for i in self.stringList()]



############################################
class SelectedFilesListView(QtWidgets.QListView):

  """
  catch mouseDoubleClickEvent to avoid edition
  """

  """not used
  def mousePressEvent(self, event):
    if verboseEvent: print "SelectedFilesListView.mousePressEvent", event.type()
    return super(SelectedFilesListView, self).mousePressEvent(event)
  """

  def mouseDoubleClickEvent(self, event):
    """avoid user entry and modify view model"""
    if verboseEvent: logger.info("SelectedFilesListView.mouseDoubleClickEvent %s" % event.type())
    #return super(SelectedFilesListView, self).mouseDoubleClickEvent(event)
    return
  

############################################
if QWebEngineView != None:

  class BrowserWebView(QWebEngineView):
    """
    As QT is an async library, you probably won't have any result 
    if you immediately try to look at the html data of your webview 
    after calling load, because it returns immediately, 
    and will trigger the loadFinished signal once the result is available. 
    You can of course try to access the html data the same way as I did in 
    the _result_available method immediately after calling load, 
    but it will return an empty page (that's the default behavior).
    """
    def __init__(self, *args):
      super(BrowserWebView, self).__init__(*args)
      #self.loadFinished.connect(self._result_available)

    def _result_available(self, ok):
      """useless TODO remove?"""
      frame = self.page().mainFrame()
      logger.info("BrowserWebView result_available")
      if verbose: logger.info("contents :\n" + str(frame.toHtml()).encode('utf-8'))

else:

  class BrowserWebView(QtWidgets.QTextBrowser):
    """
    default less functionalities browser
    if PyQt4.QtWebKit or PyQt5.QtWebEngineWidgets not installed (as salome 8.2)
    """

    def __init__(self, *args):
      super(BrowserWebView, self).__init__(*args)
      self._db = QtCore.QMimeDatabase()
 

    def load(self, url):
      typ = self._db.mimeTypeForFile(url.path()).name()
      
      if "html" in typ:
        self.setSource(url)
        return

      if "image" in typ:
        src = """\
<html>
<img src="%s" />
</html>""" % url.path()
        self.setHtml( src )
        return

      #self.mimeSourceFactory().setExtensionType("txt", "text/plain")
      #with open(url, "r") as f: contents = f.read()
      #self.setHtml( contents )
      logger.warning("QTextBrowser mimetype unknown: %s" % typ)



############################################
class ComboBoxHighlightTheme(QtWidgets.QComboBox):
  def __init__(self, *args):
    super(ComboBoxHighlightTheme, self).__init__(*args)

    if platform.system() == "windows":
      res = ["No Highlight themes for windows, are for linux"]
      self.addItems(res)
      self.setValue(res[0])
      return

    if os.path.isdir(_themesDir):
      cmd = "ls %s" % _themesDir
      rc = UXYZ.Popen(cmd)
      res_out = rc.getValue()
      res = res_out.split()
      res = [i.replace(".theme", "") for i in res if ".theme" in i]
    else:
      logger.warning("""highlight themes not found: '%s'
fix it installing package (apt-get install highlight)""" % _themesDir)
      res = ["No Highlight themes"]
    if verbose: logger.info("ComboBoxHighlightTheme themes:\n%s" % res)
    self.addItems(res)
    self.setValue(res[0])

  def setValue(self, value):
    # if verbose: print("ComboBoxHighlightTheme.setValue", str(value), type(value))
    index = self.findText(str(value))
    if index == -1:
      logger.warning("ComboBoxHighlightTheme theme not found: '%s'" % str(value))
      self.setCurrentIndex(0)
    else:
      self.setCurrentIndex(index)
    
  def getValue(self):
    """
    method implemented to avoid have to use equivalent methods with other names:
    
    - QLineEdit.text()
    - QComboBox.currentText()
    - QSpinBox/QDoubleSpinBox.value()
    - etc...
    """
    return str(self.currentText())


############################################
class ChooseHighlightTheme(QtWidgets.QDialog):
  
  highlightThemeChanged = QtCore.pyqtSignal(str) #when style changed
  
  def __init__(self, *args):
    super(ChooseHighlightTheme, self).__init__(*args)
    self.combo = ComboBoxHighlightTheme() #QtWidgets.QWidget() #
    self.prefix = "Current theme:"
    self.label = QtWidgets.QLabel(self.prefix)
    layout = QtWidgets.QVBoxLayout()
    layout.addWidget(self.label)
    layout.addWidget(self.combo)
    self.setLayout(layout)
    self.combo.activated[str].connect(self.onActivated)
    self.setWindowTitle('Choose Highlight theme')
    self.resize(300,70)
    
  def onActivated(self, text):
    self.highlightThemeChanged.emit(text)
    if verboseEvent: logger.info("emit highlightThemeChanged: '%s'" % text)
    self.close()

  def getValue(self):
    return self.combo.getValue()

  def setValue(self, value):
    return self.combo.setValue(value)


############################################
class TextView(QtWidgets.QTextEdit): #or QtWidgets.QPlainTextEdit
  """
  viewer as read only syntax highliter
  """
  def __init__(self, *args):
    super(TextView, self).__init__(*args)
    self.setContextMenuPolicy(QtCore.Qt.CustomContextMenu)
    self.customContextMenuRequested.connect(self.contextMenuSlot)
    self.currentPath = ""
    editor = UXYZ.getEditor(basename=True)
    self.actions = [
      _createAction(self, '%s edit' % editor, None, 'edit by %s' % editor, self.editBy, 'edit'), ]

    if platform.system() != "Windows":
      self.comboBoxTheme = ChooseHighlightTheme(self)
      self.comboBoxTheme.setValue("acid") #night
      self.actions.append(
        _createAction(self, 'highlight theme', None, 'change syntax highlighting theme view', self.changeThemeAction, 'highlight'))
      self.comboBoxTheme.highlightThemeChanged.connect(self.on_changeTheme)

    self.fontName = "Monospace"
    font = QtGui.QFont(self.fontName)
    font.setPointSize(9)
    #print font.family(), font.fixedPitch() #is monospaced?
    self.setFont(font)


  def contextMenuSlot(self):
    """create menu from StandardContextMenu and appending actions"""
    #warning if '_menu' becomes 'self._menu' generates seg fault when main widget close
    #seems createStandardContextMenu() owner problem
    _menu = self.createStandardContextMenu() #copy/paste/select not everytimes identical
    _menu.addSeparator()
    for action in self.actions: 
      _menu.addAction(action)
    _menu.exec_(QtGui.QCursor.pos())

  def changeThemeAction(self):
    self.comboBoxTheme.exec_() #TODO QtGui.QCursor.pos()

  def on_changeTheme(self):
    if verboseEvent: logger.info("change theme %s" % self.comboBoxTheme.getValue())
    self.openFileHighlight(self.currentPath)

  def editBy(self):
    if not os.path.isfile(self.currentPath):
      QtWidgets.QMessageBox.warning(self,"warning","Problem not a file:\n '%s'" % self.currentPath)
      return
    cmd = "%s %s &" % (UXYZ.getEditor(), self.currentPath)
    # print("begin %s '%s'" % (__file__, cmd))
    # rc = UXYZ.Popen(cmd, logger=logger) # wait end
    res = SP.Popen(cmd, shell=True)
    # print("end '%s'" % str(res))

  def openFile(self, nameFile):
    """open ascii as it without highlight"""
    if nameFile == "": return #cancel
    if isLibraryFile(nameFile):
      self.openFileLibrary(nameFile)
      return
    if not isTextFileMimetypes(nameFile):
      QtWidgets.QMessageBox.warning(self,"warning","Problem not a mimetype text file:\n%s" % nameFile)
      return
    try:
      if os.path.exists(nameFile):
        #source accent is code python with #coding=utf-8
        self.clear()
        self.setPlainText(open(nameFile, 'r').read())
        self.currentPath = str(nameFile)
        sb = self.verticalScrollBar() #scroll to the top 
        sb.triggerAction(sb.SliderToMinimum)
      else:
        self.editView.insertLine("Problem inexisting file:\n"+nameFile, "Red")
    except:
      QtWidgets.QMessageBox.warning(self,"warning","Problem reading or encoding file (as not ascii, for example):\n"+nameFile)

  def openFileLibrary(self, nameFile):
    """write contents of library without highlight"""
    if nameFile == "": return #cancel
    if not isLibraryFile(nameFile):
      QtWidgets.QMessageBox.warning(self,"warning","Problem not a library file:\n%s" % nameFile)
      return
    if True:#try:
      if os.path.exists(nameFile):
        #source accent is code python with #coding=utf-8
        cmd = """\
echo '****************** ldd -v ******************'
ldd -v %s
echo '****************** nm -gC ******************'
nm -gC %s""" % (nameFile, nameFile)
        rc = UXYZ.Popen(cmd, logger=logger)
        contents = rc.getValue()
        self.clear()
        self.setPlainText(contents)
        self.currentPath = str(nameFile)
        sb = self.verticalScrollBar() #scroll to the top 
        sb.triggerAction(sb.SliderToMinimum)
      else:
        self.editView.insertLine("Problem inexisting file:\n"+nameFile, "Red")
    else: #except: 
      QtWidgets.QMessageBox.warning(self,"warning","Problem reading file: "+nameFile)

  def openFileHighlight(self, nameFile):
    """
    try to highlight, else open ascii openFile()
    avoid highlight big file
    """
    if nameFile == "": return #cancel
    if not os.path.exists(nameFile):
      QtWidgets.QMessageBox.warning(self, "warning", "Problem inexisting file:\n  '%s'" % nameFile)
      return
    if isIradinaOutputFile(nameFile):
      return self.openFile(nameFile)
    if isLibraryFile(nameFile):
      self.openFileLibrary(nameFile)
      return      
    if not isTextFileMimetypes(nameFile):
      QtWidgets.QMessageBox.warning(self,"warning","Problem not a mimetype text file:\n  '%s'" % nameFile)
      return
    size = os.path.getsize(nameFile)
    if size > 1000000: 
      logger.warning("big file are not highlighted: %i" % size)
      return self.openFile(nameFile)
    theme = self.comboBoxTheme.getValue()
    if theme == "NoHighlightThemes": #no highlight package
      return self.openFile(nameFile)
    if True: #try:
      #source accent is code python with #coding=utf-8
      self.clear()
      # cmd = "highlight -i %s --out-format=html --style %s --font-size=8 --force" % (nameFile, theme)
      # cmd = "highlight -i %s --out-format=html --style %s --font-size=8" % (nameFile, theme)
      if isIradinaInFile(nameFile): # as windows .INI file
        cmd = "highlight -i %s --out-format=html --style %s --font-size=8 --syntax=INI" % (nameFile, theme)
      elif isIniExtensionFile(nameFile): # as windows .INI file
        cmd = "highlight -i %s --out-format=html --style %s --font-size=8 --syntax=INI" % (nameFile, theme)
      else:
        cmd = "highlight -i %s --out-format=html --style %s --font-size=8" % (nameFile, theme)
      rc = UXYZ.Popen(cmd, logger=logger)
      res_out = rc.getValue()
      if not rc.isOk():
        self.openFile(nameFile)
      else:
        htmlStr = rc.getValue()
        with open("/tmp/tmp.html", "w") as f: f.write(htmlStr)
        self.clear()
        self.insertHtml(htmlStr)
        self.currentPath = str(nameFile)
        sb = self.verticalScrollBar() #scroll to the top 
        sb.triggerAction(sb.SliderToMinimum)
    else: #except: 
      QtWidgets.QMessageBox.warning(self, "warning", "Problem reading file:\n  '%s'" % nameFile)

  def openFileHtml(self, nameFile):
    """open file as ascii html yet without highlight"""
    if nameFile == "": return #cancel
    if not os.path.exists(nameFile):
      QtWidgets.QMessageBox.warning(self, "warning", "Problem inexisting file:\n  '%s'" % nameFile)
      return
    if not isTextFileMimetypes(str(nameFile)):
      QtWidgets.QMessageBox.warning(self,"warning","Problem not a mimetype text file:\n %s" % nameFile)
      return
    try:
      #source accent is code python with #coding=utf-8
      self.clear()
      self.appendHtml(open(nameFile, 'r').read())
      sb = self.editView.verticalScrollBar() #scroll to the top 
      sb.triggerAction(sb.SliderToMinimum)
    except: 
      QtWidgets.QMessageBox.warning(self,"warning","Problem reading file:\n  '%s'" % nameFile)

def is_Highlight():
  if os.path.isdir(_themesDir):
    return True
  else:
    return False

############################################
class FileSystemModelViewerWidget(FileSystemModelWidget):
  
  """
  FileSystemModelWidget plus right third widget as quick textViewer & webViever of file
  
  pseudo modal dialog  with Cancel/Apply button use. 
  see: self._lock and self.lockForChangeModelController()
  """

  endDialog = QtCore.pyqtSignal(object)  # TODO if needed in future...
  plotFileSignal =  QtCore.pyqtSignal(object)  # may be action for plot matplolib from contents of file

  def __init__(self, *args):
    super(FileSystemModelViewerWidget, self).__init__(*args)
    self.splitter3 = QtWidgets.QSplitter(QtCore.Qt.Vertical) # textedit and/or webview
    self.editView = TextView(self)
    self.splitter3.addWidget(self.editView)
    self.webView = BrowserWebView()
    self.splitter3.addWidget(self.webView)
    self.splitter.addWidget(self.splitter3)
    self.selectedModel = SelectedFilesListModel([])
    self.selectedView = SelectedFilesListView()
    self.selectedView.setModel(self.selectedModel)
    self.is_Highlight = is_Highlight() # TODO test if not present
    #self.layout1.addWidget(self.selectedView, 1)
    #self.layout1.addWidget(self.editView.comboBoxTheme, 1)

    self.fileView.setItemDelegate(FileItemDelegate(self))
    
    self.resize(1400,800)
    self.splitter.setSizes([400, 400, 600])
    self.splitter3.setSizes([200, 0]) #hidden web view empty
    self.fileView.setColumnWidth(0, 200) #name
    self.fileView.setColumnWidth(1, 100) #size
    self.fileView.setColumnWidth(2, 1) #type
    self.fileView.setColumnWidth(3, 100) #date modified

    self.editView.setReadOnly(True)
    self.editView.setTextInteractionFlags(QtCore.Qt.TextSelectableByMouse | QtCore.Qt.TextSelectableByKeyboard)

    #self.selectedView.doubleClicked.connect(self.on_selectedView_doubleClicked) #for example for futur
    
    #context menu
    self.selectedView.setContextMenuPolicy(QtCore.Qt.CustomContextMenu)
    self.selectedView.customContextMenuRequested.connect(self.openMenuSelectedView)

    self.actionsOnSelectedView = [
      _createAction(self, "clear item", None, 'clear selected item', self.clearOneSelectedModel, 'clear'),
      _createAction(self, "clear all", None, 'clear all items', self.clearAllSelectedModel, 'clear'),
    ]
    self.menuSelectedView = QtWidgets.QMenu()
    for action in self.actionsOnSelectedView:
      self.menuSelectedView.addAction(action)

    #lock pseudo modal dialoq use of widget
    self._lock = QtCore.QMutex()
    self._lockDone = False                     #mutex protected please
    self._strForWhatDefault = "selected files"
    self._strForWhat = self._strForWhatDefault #mutex protected please
    self._methodOnApply = None                 #mutex protected please
    self._cancelApply = QtWidgets.QGroupBox()
    self._cancelApply.setTitle(self._strForWhat)
    #self._cancelApply.setFlat(True)
    self._cancel = QtWidgets.QPushButton("Cancel")
    self._cancel.setMinimumWidth(20)
    self._cancel.clicked.connect(self._onCancel)
    self._apply = QtWidgets.QPushButton("Apply")
    self._apply.setMinimumWidth(20)
    self._apply.clicked.connect(self._onApply)
    CAlayout = QtWidgets.QHBoxLayout()
    CAlayout.addWidget(self._cancel)
    CAlayout.addWidget(self._apply)
    self._cancelApply.setLayout(CAlayout)
    self._cancelApply.hide()
    
    self.layout1.addWidget(QtWidgets.QLabel("selected files"))
    self.layout1.addWidget(self.selectedView, 1)
    self.layout1.addWidget(self._cancelApply)
    self.plotWindows = []
    #self._lockForSelection() #only an example use locked

    # signal connect
    self.plotFileSignal.connect(self.plotFile) # for elementary plot

  def plotFile(self, value):
    """
    value is dict {typePlot, filePlot, ...}

    example read 4 columns iradina ascii output file in pandas
    (as csv file with typeCsv = "IradinaPlot" WITHOUT header before v1.2.0)

    >> import pandas as pd
    >> f = "/volatile/wambeke/IRADINAGUI_WORKDIR/example_UO2_250_2_4_idC/output/test.disp.sum"
    >> df = pd.read_csv(f, sep=u"\t", header=None, names="x y z value".split())
    >> print(df)
    """
    import pandaspy.pandasMainWidgetXyz as PDMW
    logger.info("plotFile\n%s" % PP.pformat(value))
    mainWin = PDMW.PandasMainDialogXyz("Iradina Plot")
    aType = value["typePlot"]
    aFile = value["filePlot"]
    mainWin.setValueFromFile(aFile, typeCsv=aType)
    mainWin.hideButtons("Apply Reset Cancel".split())
    mainWin.showButtons("Quit Help".split())
    self.on_plot_iradina(mainWin)
    return

  def on_plot_iradina(self, pandasMainDialog):
    upWid = pandasMainDialog.getUpWidgetLayout() # PandasMainWidgetXyz
    tabWid = upWid.widget
    try:
      infos = tabWid.getInfosAsDict()
      if verbose: print("plot infos:\n%s" % PP.pformat(infos))
      nbCols = infos["nbCols"]
      if nbCols == 4:
        tabWid.on_plot1_click_iradina(infos)
      if nbCols == 5:
        tabWid.on_plot2_click_iradina(infos)
    except Exception as e:
      mess = """\
There is a problem for automatic iradina plot.
You have to do it yourself.
Select columns in contents and click one plot buttons."""
      wid = MBD.getMessageBoxDialog(mess=mess, exception=e)
      wid.exec_()

    self.plotWindows.append(pandasMainDialog) # accept multiples independents plot windows
    pandasMainDialog.show()

  def _lockForSelection(self):
    self.lockForChangeModelController("select files please", self._emitEndDialog)

  def _emitEndDialog(self, selectedFiles):
    logger.info("%s _emitEndDialog for emit '%s' with:\n%s" % (self._name, self._strForWhat, selectedFiles))
    self.endDialog.emit(selectedFiles)
    return True

  def lockSetDirRootPath(self, rootPath, filters=[]):
    #mutex protected please
    ok = self._lock.tryLock()
    if not ok: #do not change, in using... for else
      res = False
    else:
      self.setDirRootPath(rootPath, filters)
      self.clearAllSelectedModel()
      res = True
      self._lock.unlock()
    return res   
    
      
  def lockForChangeModelController(self, strForWhat, methodExecOnApply):
    """
    obiously all other inherited lock... Methods begins with trylock()
    """
    #mutex protected please
    ok = self._lock.tryLock()
    if not ok:
      msg = "%s locked yet for '%s'" % (self._name, self._strForWhat)
      msgMore = "\nfix 'Apply' or 'Cancel', and try again."
      QtWidgets.QMessageBox.warning(self, "warning", msg + msgMore)
      #print msg
      return False
    self._lockDone = True
    self._strForWhat = strForWhat
    self._methodExecOnApply = methodExecOnApply
    if verboseEvent: logger.info("%s locked for %s" % (self._name, self._strForWhat))
    self._cancelApply.setTitle(self._strForWhat)
    self._cancelApply.show()
    return True
 
  def _unlockWidget(self):
    #mutex protected please
    if not self._lockDone:
      logger.error("ERROR: %s unlocked yet" % self._name)
      return False
    if verboseEvent: logger.info("%s unlocked for '%s'" % (self._name, self._strForWhat))
    self._strForWhat = self._strForWhatDefault
    self._methodExecOnApply = None
    self._cancelApply.hide()
    self._lockDone = False
    self._lock.unlock()
    return True

  def _onCancel(self):
    #mutex protected please
    if self._strForWhat == self._strForWhatDefault:
      logger.error("ERROR: %s unlocked yet on Cancel" % self._name)
      return False
    self._unlockWidget()
    return False
    
  def _onApply(self):
    #mutex protected please
    if self._strForWhat == self._strForWhatDefault:
      logger.error("ERROR: %s unlocked yet on Apply" % self._name)
      return False
    files = self.selectedModel.getStrList()
    filesEnvVar = [self.toPathWithEnvVar(i) for i in files]
    if verboseEvent: logger.info("%s 'Apply' for '%s' with\n%s" % (self._name, self._strForWhat, filesEnvVar))
    self._methodExecOnApply(filesEnvVar)
    self._unlockWidget()
    return True
   
  def resetViewStretchForWeb(self):
    s = self.splitter3.sizes()
    if verbose: logger.info("reset sizes web %s" % s)
    if s[1] < 20:
      self.splitter3.setSizes([100, 1000])

  def resetViewStretchForAscii(self):
    s = self.splitter3.sizes()
    if verbose: logger.info("reset sizes ascii %s" % s)
    if s[0] < 20:
      self.splitter3.setSizes([1000, 100])

  def on_selectedView_doubleClicked(self):
    logger.info('on_selectedView_doubleClicked')

  def openMenuSelectedView(self, position):
    #https://wiki.python.org/moin/PyQt/Creating%20a%20context%20menu%20for%20a%20tree%20view
    #print("openMenuSelectedView")
    index = self.selectedView.selectionModel().currentIndex()
    self.selectedModelIndexCurrent = index
    """
    menu = QtWidgets.QMenu()
    action = QtWidgets.QAction("clear item", self)
    action.triggered.connect(self.clearOneSelectedModel)
    menu.addAction(action)
    action = QtWidgets.QAction("clear all", self)
    action.triggered.connect(self.clearAllSelectedModel)
    menu.addAction(action)
    """
    self.menuSelectedView.exec_(self.selectedView.viewport().mapToGlobal(position))

  def clearAllSelectedModel(self):
    """reset path in self.fileModel.selectedFiles"""
    if verboseEvent: logger.info("clearAllSelectedModel")
    self.selectedModel.setStringList([])
    self.dirView.refreshFileView()

  def clearOneSelectedModel(self):
    """reset path in self.fileModel.selectedFiles"""
    if verboseEvent: logger.info("clearOneSelectedModel %i" % self.selectedView.selectionModel().currentIndex().row())
    
    """cvw
    index = self.selectedModelIndexCurrent
    model = index.model()
    if model == None: return #model None in no selection
    path = model.data(index, QtCore.Qt.EditRole)
    #self.fileModel.userSelectedFiles.remove(path)
    newList = self.selectedModel.stringList()
    print "****hello newList",newList
    newList.remove(path)
    print "****bye newList",newList
    self.selectedModel.setStringList(newList)
    """
    self.selectedModel.removeRow(self.selectedView.selectionModel().currentIndex().row())
    self.dirView.refreshFileView()
    
     
  def on_fileView_clicked(self, index):
    """edit the file clicked"""
    path = str(self.fileModel.fileInfo(index).absoluteFilePath())
    if os.path.isdir(path):
      self.setRootFile(path)
    else:
      if verboseEvent: logger.info("FileSystemModelViewerWidget.on_fileView_clicked view: '%s'" % path)
      #self.openFile(path)
      if isWebFileMimetypes(path):
        if verboseEvent: logger.info("Web %s" % mimetypes.guess_type(path))
        self.openBrowserFileHtml(path)
      if isLibraryFile(path):
        if verboseEvent: logger.info("library %s" % mimetypes.guess_type(path))
        self.openFileLibrary(path)
      if isTextFileMimetypes(path):
        if verboseEvent: logger.info("ascii %s" % mimetypes.guess_type(path))
        self.openFileHighlight(path)
      #self.openFileHtml(path)
          
  def on_fileView_doubleClicked(self, index):
    """set/reset path in self.fileModel.selectedFiles"""
    path = str(self.fileModel.fileInfo(index).absoluteFilePath())
    actual = self.selectedModel.getStrList()
    if path in actual: #remove if present yet, else append
      actual.remove(path)
    else:
      actual.append(path)
    self.selectedModel.setStringList(actual)  
    if verboseEvent: logger.info("FileSystemModelWidget.on_fileView_doubleClicked slot on:\n%s" % actual)
    self.fileView.dataChanged(index,index) #refresh
    
  """
  def event(self, event):
    if verboseEvent: print "FileSystemModelViewerWidget.event",event
    if type(event) == QtGui.QHelpEvent: 
      if verboseEvent: print "help event", event.type()
    return super(FileSystemModelWidget, self).event(event)
  """

  def openFileHighlight(self, nameFile):
    if platform.system() == "Windows" or self.is_Highlight == False:
      self.editView.openFile(nameFile) # no Highlight on windows or may be on linux
    else:
      self.editView.openFileHighlight(nameFile)

  def openFileLibrary(self, nameFile):
    self.editView.openFileLibrary(nameFile)

  def openBrowserFileHtml(self, nameFile):
    #if not self.isWebview: return  #nothing
    if nameFile == "": return #cancel
    if not os.path.exists(nameFile):
      QtWidgets.QMessageBox.warning(self, "warning", "Problem inexisting file:\n  '%s'" % nameFile)
      return
    ok = isTextFileMimetypes(str(nameFile)) or isImageFileMimetypes(str(nameFile))
    if not ok:
      QtWidgets.QMessageBox.warning(self,"warning","Problem not a mimetype file for browser:\n  '%s'" % nameFile)
      return
    try:
      logger.info("webView.load %s %s" % (type(self.webView), nameFile))
      self.webView.load(QtCore.QUrl(nameFile))
      self.nameFile = nameFile
      self.resetViewStretchForWeb()
    except:
      QtWidgets.QMessageBox.warning(self,"warning","Problem web browser reading file: "+nameFile)

def example():
  app = OnceQApplication([])
  dialog = FileSystemModelViewerWidget()
  #rootPath = os.getcwd()
  rootPath = os.path.join("$URANIESYS", "share", "uranie")
  dialog.setDirRootPath(rootPath, filters=[])
  dialog.show()
  app.exec_()

if __name__ == '__main__':
  example()
  pass

