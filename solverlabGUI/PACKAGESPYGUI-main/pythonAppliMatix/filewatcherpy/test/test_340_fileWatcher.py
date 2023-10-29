#!/usr/bin/env python
# -*- coding: utf-8 -*-

# %% LICENSE_SALOME_CEA_BEGIN
# see PACKAGESPY/LICENCE file
# %% LICENSE_END


import os
#import glob
import unittest
import filewatcherpy.fileWatcher as FWR
from salomepy.onceQApplication import OnceQApplication

verbose = False

homeDir = os.path.realpath(os.getenv("HOME"))

testDir = os.path.split( os.path.realpath(__file__) )[0]

class TestCase(unittest.TestCase):

  def test_010(self):
    #as bash interpretation for "echo theParameterOfgetRealPath"
    currenDir = os.getcwd()
    self.assertEqual( FWR.getRealPath(""), currenDir)
    self.assertEqual( FWR.getRealPath("."), currenDir)
    self.assertEqual( FWR.getRealPath("$HOME"), homeDir)
    self.assertEqual( FWR.getRealPath("$HOME/toto"), homeDir + "/toto")
    self.assertEqual( FWR.getRealPath("$HOMEtoto"), currenDir + "/$HOMEtoto") #inexisting env var ${HOMEtoto}
    self.assertEqual( FWR.getRealPath("${HOME}toto"), homeDir + "toto") #best interpretation


  def test_020(self):
    app = OnceQApplication()
    afw = FWR.getFileWatcher()
    aPath = "$HOME/.bashrc"
    afw.addPath(aPath)
    if verbose: afw.addPath("$HOME/inexistingggg") #generate warning 'inexisting file'
    self.assertEqual( len(afw.getFilesWatched()), 1)
    self.assertTrue( ".bashrc" in afw.getFilesWatched()[0] )
    self.assertTrue( homeDir in afw.getFilesWatched()[0] )
    self.assertTrue( aPath in afw.getOriginalPathsWatched() )
    self.assertFalse( "toto" in afw.getOriginalPathsWatched() )

    afw.removePath(aPath)
    self.assertEqual( len(afw.getFilesWatched()), 0)
    self.assertEqual( len(afw.getOriginalPathsWatched()), 0)

    afw.addPath(aPath)
    self.assertEqual( len(afw.getFilesWatched()), 1)
    self.assertTrue( ".bashrc" in afw.getFilesWatched()[0] )
    self.assertTrue( homeDir in afw.getFilesWatched()[0] )

    if verbose: afw.addPath(aPath) #generate warning 'existing yet watched file'

    if verbose: afw.removePath("toto")

    afw.removeAllPaths() #stay aPath to remove
    self.assertEqual( len(afw.getFilesWatched()), 0)
    self.assertEqual( len(afw.getDirectoriesWatched()), 0)
    self.assertEqual( len(afw.getOriginalPathsWatched()), 0)

    #if verbose: afw.deleteOnTimer(0)

  def test_030(self):
    app = OnceQApplication()
    afw = FWR.getFileWatcher()
    aPath = "$HOME"
    afw.addPath(aPath)
    self.assertEqual( len(afw.getDirectoriesWatched()), 1)
    self.assertTrue( homeDir in afw.getDirectoriesWatched()[0] )
    self.assertTrue( aPath in afw.getOriginalPathsWatched() )
    self.assertFalse( "toto" in afw.getOriginalPathsWatched() )

    afw.removePath(aPath)
    self.assertEqual( len(afw.getFilesWatched()), 0)
    self.assertEqual( len(afw.getDirectoriesWatched()), 0)
    self.assertEqual( len(afw.getOriginalPathsWatched()), 0)

    afw.addPath(aPath)
    self.assertEqual( len(afw.getFilesWatched()), 0)
    self.assertEqual( len(afw.getDirectoriesWatched()), 1)

    if verbose: afw.addPath(aPath) #generate warning 'existing yet watched directory'
    afw.removeAllPaths() #stay aPath to remove
    self.assertEqual( len(afw.getFilesWatched()), 0)
    self.assertEqual( len(afw.getDirectoriesWatched()), 0)
    self.assertEqual( len(afw.getOriginalPathsWatched()), 0)

  def test_999(self):
    #only to test manually Event :
    # fileWatcherTest.py
    # touch ~/tototmp
    if verbose:
      print("""
fileWatcherTest.py will not quit to show event in ~ directory
  so you have to kill process unittest yourself
  try for example -> 'touch ~/tototmp'
""")
      app = OnceQApplication()
      afw = FWR.getFileWatcher()
      afw.verboseEvent = True
      aPath = "$HOME"
      afw.addPath(aPath)
      app.exec_() #do not quit so have to kill process


if __name__ == '__main__':
  #verbose = True
  FWR.verbose = verbose
  unittest.main()
  pass
