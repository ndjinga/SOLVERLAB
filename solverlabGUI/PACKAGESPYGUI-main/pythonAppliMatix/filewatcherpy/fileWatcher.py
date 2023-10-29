#!/usr/bin/env python
# -*- coding: utf-8 -*-

# %% LICENSE_SALOME_CEA_BEGIN
# see PACKAGESPY/LICENCE file
# %% LICENSE_END


import os
from PyQt5 import QtCore 
import subprocess as SP

import xyzpy.loggingXyz as LOG

"""
see http://doc.qt.io/qt-4.8/qfilesystemwatcher.html
Note: The act of monitoring files and directories for modifications consumes system resources. This implies there is a limit to the number of files and directories your process can monitor simultaneously

so user have to avoid multiple instanciation of FileWatcher, use getFileWatcher() singleton
"""

logger = LOG.getLogger()
verbose = False


_theFileWatcher = [] #the singleton: avoid multiples instances: only one loop timer


################################################
class FileWatcher(QtCore.QFileSystemWatcher):
  
  """
  class for manage changes in files or directories from QtCore.QFileSystemWatcher
  theorically user have to only use singleton instance from getFileWatcher()
  append funtionnality 'interpret path with env variable as $HOME/toto etc...' 
  """  
  
  def __init__(self, *args, **kwargs):
    super(FileWatcher, self).__init__(*args, **kwargs)
    self.fileChanged.connect(self.aFileChanged)
    self.directoryChanged.connect(self.aDirectoryChanged)
    #memorize original files names with env var as "$HOME/toto"
    self._originalPaths = {} #memorize original dir/files names with env var as "$HOME/toto"
    self.verboseEvent = False

  def addPath(self, aPathWithEnvVar):
    aRealPath = getRealPath(aPathWithEnvVar)
    if os.path.isfile( aRealPath ):
      if aRealPath in self.files(): #key exists
        logger.warning("existing yet watched file: '%s'"% aRealPath)
        return
      self._originalPaths[aRealPath] = aPathWithEnvVar
      super(FileWatcher, self).addPath(aRealPath)
    elif os.path.isdir( aRealPath ):
      if aRealPath in self.directories(): #key exists
        logger.warning("existing yet watched directory: '%s'"% aRealPath)
        return
      self._originalPaths[aRealPath] = aPathWithEnvVar
      super(FileWatcher, self).addPath(aRealPath)
    else:
      logger.warning("inexisting file or directory: '%s'"% aRealPath)

  def addPaths(self, aListPathWithEnvVar):
    for i in aListPathWithEnvVar:
      self.addPath(i)

  def removePath(self, aPathWithEnvVar):
    aRealPath = getRealPath(aPathWithEnvVar)
    if os.path.isfile( aRealPath ):
      if aRealPath not in self.files(): #key not exists
        logger.warning("not currently watched file: '%s'"% aRealPath)
        return    
    elif os.path.isdir( aRealPath ):
      if aRealPath not in self.directories(): #key not exists
        logger.warning("not currently watched directory: '%s'"% aRealPath)
        return
    else:
      logger.warning("inexisting file or directory: '%s'"% aRealPath)
    if aRealPath in self._originalPaths: #in case of desynchronize  
      del self._originalPaths[aRealPath]
    super(FileWatcher, self).removePath(aRealPath)
    
  def removePaths(self, aListPathWithEnvVar):
    for i in aListPathWithEnvVar:
      self.removePath(i)

  def removeAllPaths(self):
    self._originalPaths = {}
    for i in self.files(): super(FileWatcher, self).removePath(i)
    for i in self.directories(): super(FileWatcher, self).removePath(i)
      
  def aFileChanged(self, aFile):
    """to print event if verboseEvent"""
    #print "files  watched:", self.getFilesWatched()
    if self.verboseEvent: print("FileWatcher.FileChanged: '%s'" % aFile)

  def aDirectoryChanged(self, aDir):
    """to print event if verboseEvent"""
    #print "directories  watched:", self.getDirectoriesWatched()
    if self.verboseEvent: print("FileWatcher.DirectoryChanged: '%s'" % aDir)

  def getFilesWatched(self):
    return [str(i) for i in self.files()]
       
  def getDirectoriesWatched(self):
    return [str(i) for i in self.directories()]
       
  def getOriginalPathsWatched(self):
    return sorted([str(i) for k, i in list(self._originalPaths.items())])
       
  """
  def event(self, event):
    print "event",event
    return super(FileWatcher, self).event(event)
  """
    
  def deleteOnTimer(self, deltaTime):
    if id(self) == id(_theFileWatcher[0]):
      logger.warning("dangerous on singleton _theFileWatcher: delete not done")
      return
    timer = QtCore.QTimer(self);
    timer.timeout.connect(self.deleteLater)
    timer.start(deltaTime)


################################################
def getRealPath_obsolete(aPathWithEnvVar):
  """
  resolve env variable as $HOME/toto etc...
  with subprocess shell interpretation of env var
  """
  cmd = "echo %s" % aPathWithEnvVar
  proc = SP.Popen(cmd, shell=True, stdout=SP.PIPE, stderr=SP.PIPE)
  stdout, stderr = proc.communicate()
  stdout = stdout.rstrip()
  stderr = stderr.rstrip()
  # if verbose: print("stdout: '%s'\nstderr: '%s'" % (stdout, stderr))
  try:
    res = os.path.realpath(stdout)
  except:
    raise Exception("Problem shell interpretation of path: '%s'" % aPathWithEnvVar)

  if verbose: print("getRealPath: '%s'->'%s'" % (aPathWithEnvVar, res))
  return str(res)


################################################
def getRealPath(aPathWithEnvVar):
  """
  resolve env variable as $HOME/toto etc...
  with expandvars
  """
  res = os.path.realpath(os.path.expandvars(aPathWithEnvVar))
  return res


################################################
def getFileWatcher():
  """
  use it as singleton
  """
  if len(_theFileWatcher) == 0:
    _theFileWatcher.append(FileWatcher())
  return _theFileWatcher[0]


################################################
def exampleLaunchStandalone():
  from salomepy.onceQApplication import OnceQApplication
  app = OnceQApplication()
  aFileWatcher = FileWatcher()
  aFileWatcher.verboseEvent = True
  #for i in dir(aFileWatcher): print i
  aPath = "essai0.tmp"
  aFileWatcher.addPath(aPath)
  aPaths = ["essai1.tmp", "essai2.tmp"]
  aFileWatcher.addPaths(aPaths) 
  #aFileWatcher.deleteOnTimer(2000)  
  #for i in dir(app): print i
  timer = QtCore.QTimer(app);
  timer.timeout.connect(app.quit)
  timer.start(30000)
  
  print("\n!!!! begin exampleLaunchStandalone !!!!")
  print("files to modify: ", aFileWatcher.getFilesWatched())
  print("example to modify: touch ", aFileWatcher.getFilesWatched()[0])
  app.exec_()
  print("!!!! end on timer quit !!!!")

 
################################################
if __name__ == '__main__':
  exampleLaunchStandalone()

