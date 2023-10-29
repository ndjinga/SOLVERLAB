#!/usr/bin/env python
# -*- coding: utf-8 -*-

# %% LICENSE_SALOME_CEA_BEGIN
# see PACKAGESPY/LICENCE file
# %% LICENSE_END

"""
#to lanch test standalone:
matix shell
AllTestLauncher.py $PACKAGESPY_ROOT_DIR utilsSubprocessTest.py
"""

import os
import unittest
import time

import salomepy.utilsWorkdir as UTW
import salomepy.utilsSubprocess as UTS
import xyzpy.loggingXyz as LOG

logger = LOG.getLogger()


verbose = False
myDir, myName =  os.path.split(__file__)

runningDir = UTW.getWorkdirDefault("TESTS")
if os.path.isdir(runningDir):
  if verbose: logger.info("runningDir %s existing" % runningDir)
else:
  logger.warning("runningDir %s NOT existing, create it." % runningDir)
  os.makedirs(runningDir)

  
_dgibis = {

"small_test" : """\
*23456789123456789123456789123456789123456789123456789123456789123456789
*=======================================================================
*
*       A SMALL TEST
*
*=======================================================================
*
*
* DEBUT DU FICHIER MAILLE1.DGIBI
*
TITR ‘PREMIER MAILLAGE’ ;
OPTI DIME 3 ELEM CUB8 ;
*
* POINTS
*
P1 = 0. 0. 0. ;
P2 = 10. 0. 0. ;
P3 = 10. 10. 0. ;
P4 = 0. 10. 0. ;
*
* LIGNES
*
P1P2 = P1 DROI 10 P2 ;
P2P3 = P2 DROI 10 P3 ;
P3P4 = P3 DROI 10 P4 ;
P4P1 = P4 DROI 10 P1 ;
*
* SURFACE
*
SURF1 = DALL P1P2 P2P3 P3P4 P4P1 ;
*
* VOLUME
*
VOL1 = SURF1 VOLU 10 TRANS (0. 0. 10.) ;
*
*TRACAGE
*
*TRAC SURF1 QUAL ;
*TRAC CACH VOL1 ;
*
*FIN DU FICHIER
*
OPTI SAUV FORMAT 'small_test.sauv' ;
SAUV VOL1 ;
FIN ;
""",

"small_test_bug" : """\
*23456789123456789123456789123456789123456789123456789123456789123456789
*=======================================================================
*
*       A SMALL TEST WITH BUG
*
*=======================================================================
*
*
OPTI DIME 3 ELEM CUB8 ;
*
opti rest 'inexistingFileForBug.sauv';
MESS 'I NEVER COME HERE' ;
FIN ;
""",

"small_test_infini" : """\
*23456789123456789123456789123456789123456789123456789123456789123456789
*=======================================================================
*
*       A SMALL TEST WITH LONG EXECUTION TIME
*
*=======================================================================
*
*
TITR ‘PREMIER MAILLAGE’ ;
OPTI DIME 3 ELEM CUB8 ;
*
I = 0 ;
* 1000000 is more than 10 seconds?,... depends of machine
REPETER NOMBOUCL 1000000 ;
  i = i + 1 ;
FIN NOMBOUCL ;
MESS 'RESULT I =' i ;
FIN ;
""",

}

def prepareSrcDir(name, runningDir):
  srcDir = os.path.join(runningDir, "srcDir_"  + name)
  destDir = os.path.join(runningDir, "runDir_" + name)
  UTW.rmTree(srcDir)
  UTW.makeDir(srcDir)
  #if verbose: print "srcDir", srcDir
  nameFile = os.path.join(srcDir, name + ".dgibi")
  with open(nameFile, 'w') as f:
    f.write(_dgibis[name])
  return (srcDir, destDir)

def isCastemInEnv():
  for k in list(os.environ.keys()):
    if 'castem' in k.lower(): return True
  return False
  
  
class TestCase(unittest.TestCase):
  def test_000(self):
    # small_test launch test
    if not isCastemInEnv():
      logger.warning("CASTEM environment variable not found, skip tests")
      return

    if verbose: logger.info("runningDir", runningDir)
    name = "small_test"
    srcDir, destDir = prepareSrcDir(name, runningDir)
    timeOut = 10
    resultFile = name + ".sauv"
    ok = UTS.runCastem(name, srcDir, destDir, timeOut, resultFile=resultFile)
    self.assertEqual(ok, "ok")

  def test_001(self):
    # small_test_bug launch test
    if not isCastemInEnv(): return
    name = "small_test_bug"
    srcDir, destDir = prepareSrcDir(name, runningDir)
    timeOut = 10
    resultFile = name + ".sauv"
    ok = UTS.runCastem(name, srcDir, destDir, timeOut, resultFile=resultFile)
    self.assertEqual("ko" in ok, True)
    self.assertEqual("Result file not found" in ok, True)

  def test_002(self):
    # small_test timeout launch test
    if not isCastemInEnv(): return
    name = "small_test_infini"
    srcDir, destDir = prepareSrcDir(name, runningDir)
    destDir += "_002"
    timeOut = 1
    resultFile = name + ".sauv"
    ok = UTS.runCastem(name, srcDir, destDir, timeOut, resultFile=resultFile)
    self.assertEqual("ko" in ok, True)
    self.assertEqual("timeout reached" in ok, True)

  def test_003(self):
    if not isCastemInEnv(): return
    name = "small_test_infini"
    srcDir, destDir = prepareSrcDir(name, runningDir)
    destDir += "_003" #avoid working on previous _002 destDir
    timeOut = 1
    resultFile = name + ".sauv"
    ok = UTS.runCastem(name, srcDir, destDir, timeOut, resultFile=None)
    self.assertEqual("ko" in ok, True)
    self.assertEqual("timeout reached" in ok, True)

  def test_004(self):
    if not isCastemInEnv(): return
    name = "small_test_infini"
    srcDir, destDir = prepareSrcDir(name, runningDir)
    destDir += "_004" #avoid working on previous _003 destDir
    timeOut = 100
    resultFile = name + ".sauv"
    ok = UTS.runCastem(name, srcDir, destDir, timeOut, resultFile=None)
    self.assertEqual("ok", ok)
    
  def test_999(self):
    UTW.rmTree(runningDir)


if __name__ == '__main__':
  verbose = True   #verbose if human eye
  unittest.main()

