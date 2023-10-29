#!/usr/bin/env python
# -*- coding: utf-8 -*-

# %% LICENSE_SALOME_CEA_BEGIN
# see PACKAGESPY/LICENCE file
# %% LICENSE_END

"""
goal is windows where echoes colored trace messages of subprocess execution
"""

"""
example of subprocess.Popen

>>> import subprocess as SP
>>> proc = SP.Popen(["python", "-u", "-c", "print 'hello word from python Popen!'"],
>>>                 stdin=SP.PIPE, stdout=SP.PIPE, stderr=SP.PIPE)
>>> proc = SP.Popen("ls -alt", shell=True, stdin=SP.PIPE, stdout=SP.PIPE, stderr=SP.PIPE)

"""

import os
import sys
import subprocess as SP
import pprint as PP

from PyQt5 import QtGui, QtCore, QtWidgets

from datetime import datetime
from salomepy.threadWorkerForWidgetEdit import ThreadWorkerForWidgetEdit
import salomepy.iconsUser as IUSR
import xyzpy.loggingXyz as LOG

logger = LOG.getLogger()

verbose = False
verboseEvent = verbose

_do_only_one = []

class QTextEditForLog(QtWidgets.QPlainTextEdit):
  """
  implements a widget QPlainTextEdit, 
  with 2 workers as 2 threading processes
  waiting for stdout and stderr of a subprocess.
  colored echoes in window
  """
  
  # Define signals
  trigger = QtCore.pyqtSignal(str, str)
  finishPopen = QtCore.pyqtSignal(str) #when finish in popen
  index = [0] #need unambigous name...
  
  def __init__(self):
    super(QTextEditForLog, self).__init__()
    self.nameFile = None #a command popen log is not an opened file
    self.setObjectName("QTextEditForLog_%i" % self.index[0])
    self.setWindowTitle(self.objectName())
    self.index[0] += 1
    self.saveFileExt = ".log"

    """
    font = QtGui.QFont()
    font.setFamily("Courier")
    font.setBold(False)
    font.setWeight(100)
    self.setFont(font)
    font.setStyleHint(font.Monospace);
    self.setStyleSheet('font: 9pt "%s";' % self.fontName )
    """

    """
    font = QtWidgets.QFontDialog.getFont()
    try:  # return tuple sometimes (<PyQt5.QtGui.QFont>, True)
      font = font[0]
    except:
      pass
    logger.info("font family %s fixed pitch %s" % (font.family(), font.fixedPitch()))
    """

    # self.fontName = "Century Schoolbook L"   # proportional
    self.fontName = "Monospace"
    # self.fontName = "Courier"
    # self.fontName = "FreeMono"
    self.fontColor = "black"
    self.fontSize = 9
    self.fontStyle = "normal"
    # font = QtGui.QFont(self.fontName)
    QFD = QtGui.QFontDatabase
    font = QFD.systemFont(QFD.FixedFont)

    font.setPointSize(self.fontSize)
    # https://stackoverflow.com/questions/10977259/qt-fixed-width-font
    font.setFixedPitch(True) # do not work if Century Schoolbook L
    # logger.info("font family %s fixed pitch %s" % (font.family(), font.fixedPitch()))
    self.setFont(font)

    self.setLineWrapMode(self.NoWrap)
    self.resize(600,500)
    self.trigger.connect(self.insertLine)
    self._Workers = {}
    
    self.workDir = None #if popen log
    self.currentDir = "" #if open file
    self.initializeWorkdir()
    
    self.popenCommand = None
    self.withEditor = os.getenv('EDITOR4LOG')
    if self.withEditor == None:
      try:
        import iradinapy.iradinaSettings as ISET
        edit = ISET.getVar("editor")
      except:
        edit = "pluma"
        msg = "default editor is %s" % edit
        if msg not in _do_only_one:
          logger.debug(msg)
          _do_only_one.append(msg)

      self.withEditor = edit

    clearAction=QtWidgets.QAction("Clear All", self)
    clearAction.triggered.connect(self.clear)

    editAction=QtWidgets.QAction("Edit with "+self.withEditor, self)
    editAction.setIcon(IUSR.getIconFromName("editor"))
    editAction.triggered.connect(self.bestEdit)

    saveAction = QtWidgets.QAction('Save', self)
    saveAction.setShortcut('Ctrl+S')
    saveAction.setStatusTip('Save current file')
    saveAction.triggered.connect(self.saveFile)
    
    openAction = QtWidgets.QAction('Open', self)
    openAction.setShortcut('Ctrl+O')
    openAction.setStatusTip('Open a file')
    openAction.triggered.connect(self.openKnownFile)
    
    editContextMenu = self.createStandardContextMenu()
    
    editContextMenu.addAction(clearAction)
    #self.editContextMenu.addSeparator()
    action0=editContextMenu.actions()[0]
    editContextMenu.insertAction(action0, openAction)
    editContextMenu.insertAction(action0, saveAction)
    editContextMenu.insertAction(action0, editAction)
    editContextMenu.insertSeparator(action0)
    
    for a in editContextMenu.actions(): a.setEnabled(True)
    
    self.editContextMenu = editContextMenu
    #Erreur de segmentation: need setparent
    editContextMenu.setParent(self)
 
    pal=self.palette()
    pal.setColor(pal.Base, QtGui.QColor(230,230,230))
    pal.setColor(pal.Text, QtGui.QColor(0,0,0))
    
    self.setPalette(pal)
    self.finishPopen.connect(self.finishPopenAction)
    
    self.signal2Emit = None
    pass

  def initializeWorkdir(self):
    """
    initialize user logs directory $WORKDIR4LOG if not existing
    """
    workDir=os.getenv("WORKDIR4LOG")
    if workDir == None: 
      workDir=os.path.join("/tmp",  os.getenv('USERNAME'))
    workDir = os.path.realpath(workDir)
    os.environ["WORKDIR4LOG"] = workDir
    if not os.path.exists(workDir):
      logger.warning('create inexisting $WORKDIR4LOG ' + workDir)
      os.makedirs(workDir)
    logger.debug("QTextEditForLog logs files in $WORKDIR4LOG '%s'" % workDir)
    self.workDir = workDir

  def closeEvent(self, event):
    if verboseEvent:
      print(self.objectName()+".closeEvent with popen in thread", len(list(self._Workers.keys())))
    for _, (iout, ierr) in list(self._Workers.items()): 
      if iout != None: iout.stop() #reentrance on thread event
      if ierr != None: ierr.stop() #reentrance on thread event
    self._Workers = {} #allow garbage collection of threads
    return super(QTextEditForLog, self).closeEvent(event)
    
  def fromUtf8(self,  value):
    """useless pyqt5 returns unicode, QString obsolete"""
    #print type(value), value
    if sys.version_info[0] < 3:
      return u"%s" % value
    else:
      return str(value)
    
  def getCurrentFontFamilies(self):
    res = []
    fontDatabase = QtGui.QFontDatabase()
    for family in fontDatabase.families():
      if verbose: print(self.fromUtf8(family))
      res.append(self.fromUtf8(family))
    return res

  def contextMenuEvent(self, event):
    if verboseEvent:
      print("DEBUG: " + self.objectName() + " contextMenuEvent")
    self.editContextMenu.exec_(event.globalPos())

  """
  Usage: kate [Qt-options] [KDE-options] [KDE-tempfile-options] [options] [URL] 
  Kate - Advanced Text Editor
  Generic options:
    --help                    Show help about options
    --help-qt                 Show Qt specific options
    --help-kde                Show KDE specific options
    --help-kde-tempfile       Show KDE-tempfile specific options
    --help-all                Show all options
    --author                  Show author information
    -v, --version             Show version information
    --license                 Show license information
    --                        End of options
  Options:
    -s, --start <name>        Start Kate with a given session
    -u, --use                 Use a already running kate instance (if possible)
    -p, --pid <pid>           Only try to reuse kate instance with this pid
    -e, --encoding <name>     Set encoding for the file to open
    -l, --line <line>         Navigate to this line
    -c, --column <column>     Navigate to this column
    -i, --stdin               Read the contents of stdin
  Arguments:
    URL                       Document to open
  """
  def bestEdit(self):
    if self.nameFile==None:
      self.nameFile=self.getDefaultNameFile()
      self.saveFile()
    #2> /dev/null because sometimes warnings editors: Can't load fallback CSS etc...
    cmd = self.withEditor + " " + self.nameFile + " 2> /dev/null &"
    if self.withEditor=="gedit":
      cmd = self.withEditor + " " + self.nameFile + " 2> /dev/null &"
    if self.withEditor=="kate":
      # do not work on centos6.4 cmd = self.withEditor + " --start cassis --use " + self.nameFile + " &"
      cmd = self.withEditor + " " + self.nameFile + " 2> /dev/null &"
    proc = SP.Popen(str(cmd), shell=True)

  def saveFile(self):
    if self.nameFile==None:
      nameFile = QtWidgets.QFileDialog.getSaveFileName(self, 'Save File', self.currentDir)[0]
    else:
      nameFile = self.nameFile
    if nameFile=="": return #cancel
    realPath = os.path.realpath(nameFile)
    dirName = os.path.dirname(realPath)
    if not os.path.exists(dirName): os.makedirs(dirName)
    filedata = self.document().toPlainText()
    # logger.debug("saveFile %s type %s" % (realPath, type(filedata)))
    with open(realPath, 'w') as f:
      if sys.version_info[0] < 3:
        f.write(filedata.encode('utf8')) # skip colors formatting
      else:
        f.write(filedata)  # TODO skip colors formatting
    self.nameFile = nameFile

  def getDefaultNameFile(self):
    ext = datetime.now().strftime("%Y_%m_%d-%H_%M_%S")
    return os.path.join(self.workDir, "log_" + ext + self.saveFileExt)

  def setNameFile(self, nameFile):
    """needs nameFile with path"""
    self.nameFile = os.path.realpath(nameFile)
 
  def setCurrentDir(self, currentDir):
    """needs arg currentDir for path of log file"""
    self.currentDir = currentDir
    logger.debug("currentDir of open/save file " + self.currentDir)
 
  def openKnownFile(self):
    if self.nameFile==None:
      nameFile = QtWidgets.QFileDialog.getOpenFileName(self, 'Open File', self.currentDir)[0]
      if nameFile=="": return #cancel
    else:
      nameFile = self.nameFile
    self.openFile(nameFile)

  def openFile(self, nameFileIni):
    nameFile = os.path.expandvars(nameFileIni)
    if nameFile == "": return #cancel
    if '.png' in nameFile:
      QtWidgets.QMessageBox.warning(self,"warning","Problem not a text file: "+nameFile)
      return
    try:
      if os.path.exists(nameFile):
        #source accent is code python with #coding=utf-8
        self.clear()
        self.setPlainText(self.fromUtf8(open(nameFile, 'r').read()))
        self.nameFile = nameFile
        self.setCurrentDir(os.path.dirname(os.path.realpath(nameFile)))
      else:
        #self.insertLine("Problem inexisting file:\n"+nameFile, "Red")
        QtWidgets.QMessageBox.warning(self,"warning","Problem inexisting file:\n'%s'" % nameFile)
    except: 
      QtWidgets.QMessageBox.warning(self,"warning","Problem reading file: "+nameFile)
  
  def moveCursor(self, operation, mode=QtGui.QTextCursor.MoveAnchor):
     """
     TODO do not work...
      move the cursor. operation could be:
      
      - QtGui.QTextCursor.End
      - QtGui.QTextCursor.Left
      - QtGui.QTextCursor.Right
      - QtGui.QTextCursor.EndOfLine
      - QtGui.QTextCursor.LineUnderCursor
     """
     logger.warning("TODO moveCursor do not work, have to debug")
     cursor = self.textCursor()
     cursor.movePosition(operation, mode)
     self.setTextCursor(cursor)

  def insertLine(self, line, color=None, fontName=None, fontSize=None, fontStyle=None):
     """insert one standard line without forbidden tag xml <>"""
     
     """
     see http://en.wikipedia.org/wiki/Font_family_%28HTML%29
     #example of format:
     <pre style="font-family: times, serif; font-size:14pt; font-style:italic; color:red">
     Sample text formatted with inline CSS.
     </pre>
     warning <p> strip multiples whitespaces
             <pre> do not strip multiples whitespaces
     """
     
     if line[-1] == "\n": line = line[:-1] #appendHtml do one lf
     
     col, nam, siz, sty = (self.fontColor, self.fontName, self.fontSize, self.fontStyle)
     if color != None: col = self.fromUtf8(color)
     if fontName != None: nam = self.fromUtf8(fontName)
     if fontSize != None: siz = self.fromUtf8(fontSize)
     if fontStyle != None: sty = self.fromUtf8(fontStyle)
     
     #http://www-sul.stanford.edu/tools/tutorials/html2.0/whitespace.html
     aTag='<pre style="font-family: %s; font-size:%spt; font-style:%s; color:%s">%s</p>'
     #source accent is code python with #coding=utf-8
     a=aTag % (nam, siz, sty,  color, self.fromUtf8(line) )
     self.appendHtml( a )
     self.ensureCursorVisible()
     
  def insertText(self, aText):
     """
     insert text as it
     """
     self.appendPlainText( aText )
     self.ensureCursorVisible()
     
  def insertTextColor(self, aText, color="Black"):
     """could be risky if text have <xml tags> expressions"""
     if "<" in aText:
       self.insertText(aText)
     else:
       self.insertLine(aText, color)

  def indent(self, aStr):
    """simple indentate all lines"""
    ind = "  " # indentation
    s = str(aStr) # if byte
    s = s.split("\n")
    s = ("\n" + ind).join(s)
    return ind + s

  def createAndStartWorkers(self, proc):
    """
    Workers are 2 threading processes.
    Created to wait for stdout and stderr of subprocess, 
    and colored echoes in self (wich is widget QPlainTextEdit).
    
    proc is processus, the result of a call subprocess.Popen()
    """
    stdoutWorker = ThreadWorkerForWidgetEdit(proc.stdout, self, "Black")
    stderrWorker = ThreadWorkerForWidgetEdit(proc.stderr, self, "Red", emit=False)
    # warning 2017 : I never seen finishPopen on stderr don't know why...
    # so set stdout stderr both for garbage collecting on stdoutWorker signal
    self._Workers[stdoutWorker.objectName()] = (stdoutWorker, stderrWorker)
    if verbose: 
      print("%s._Workers %s" % (self.objectName(), PP.pformat(self._Workers)))
    stdoutWorker.start()
    stderrWorker.start()
    
  def launchIntoPopen(self, cmd, signal=None, cwd=None):
    """launch as background"""
    self.nameFile = None #a command is not a file
    self.insertText(cmd)
    if type(cmd)==list:
      proc = SP.Popen(cmd, shell=False, stdout=SP.PIPE, stderr=SP.PIPE, cwd=cwd)
    else:
      proc = SP.Popen(cmd, shell=True, stdout=SP.PIPE, stderr=SP.PIPE, cwd=cwd)
    self.popenCommand = cmd
    self.signal2Emit = signal
    if verbose: print("createAndStartWorkers cmd:\n%s" % self.indent(cmd))
    self.createAndStartWorkers(proc)
    
  def finishPopenAction(self, threadName):
    if self.signal2Emit != None:
      self.signal2Emit.emit()
    
    #if "castem" in self.popenCommand: print "TODO refresh", self.parent()
    if verbose:
      print("finishPopenAction on %s %s" % (self.objectName(), threadName))
    
    """save log file only at end popen"""
    #if name != self.objectName(): print "!!!!problem unexpected name from thread:", name, self.objectName()
    
    #allow garbage collection ?of threads?
    try:
      #print "-----del", str(threadName)
      #del self._Workers[str(threadName)] #reentrance, cause problem?
      self._Workers[str(threadName)] = (None, None) #better
    except:
      #print "\nProblem del thread\n", str(threadName), sorted(self._Workers.keys())
      pass

    if self.popenCommand != None: #2 events: from stdout, and from stderr
      self.popenCommand = None
      self.nameFile = self.getDefaultNameFile()
      self.saveFile()

if __name__ == "__main__":
  app = QtWidgets.QApplication([])
  edit = QTextEditForLog()
  edit.show()
  for i in edit.getCurrentFontFamilies():
    #first arrived, first shown
    edit.insertLine(i, fontName=i)
  app.exec_()

