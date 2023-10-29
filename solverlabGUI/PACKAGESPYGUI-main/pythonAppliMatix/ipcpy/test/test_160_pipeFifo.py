#!/usr/bin/env python
# -*- coding: utf-8 -*-

# %% LICENSE_SALOME_CEA_BEGIN
# see PACKAGESPY/LICENCE file
# %% LICENSE_END

"""
Code during executing thread
From the Open Group page for the open() function
O_NONBLOCK
  When opening a FIFO with O_RDONLY or O_WRONLY set:
  If O_NONBLOCK is set:
      An open() for reading only will return without delay. An open()
      for writing only will return an error if no process currently
      has the file open for reading.

  If O_NONBLOCK is clear:
      An open() for reading only will block the calling thread until a
      thread opens the file for writing. An open() for writing only
      will block the calling thread until a thread opens the file for
      reading.
"""

import os
import stat
import time
import pprint as PP
import subprocess as SP
import unittest

import ipcpy.pipeFifo as PFIFO
import xyzpy.loggingXyz as LOG

logger = LOG.getLogger()

_aDir = os.path.realpath(os.path.join("/tmp", os.getenv("USER"), "testpipefifo"))
_aPipeName = os.path.realpath(os.path.join(_aDir, "pipefifo_test_001"))

verbose = False

class TestCase(unittest.TestCase):

  def removePipe(self):
    try:
      os.remove(_aPipeName)
    except:
      pass

  def test_000(self):
    # done firstly
    self.removePipe()
    self.assertFalse(os.path.exists(_aPipeName))
    if verbose:
      LOG.pushLevel("INFO")
      logger.warning("logger.setLevel INFO")

  def test_999(self):
    # done lastly
    if verbose:
      LOG.popLevel("INFO")
      logger.warning("logger.setLevel as previous")

  def test_010(self):
    if verbose: logger.info("test_010")
    a = PFIFO.PipeFifo()
    self.assertEqual(a.getNameFifo(), None)
    self.assertFalse(a.isPipeFifo())
    a.setNameFifo(_aPipeName)
    self.assertTrue(_aPipeName in str(a))
    self.assertTrue(a.isPipeFifo())
    self.assertFalse(PFIFO.isFilePipeFifo(_aDir))
    self.assertFalse(PFIFO.isFilePipeFifo(__file__))

  def test_020(self):
    if verbose: logger.info("test_020")
    a = PFIFO.PipeFifo(_aPipeName)
    self.assertTrue(a._testPipeFifoBlocking())
    self.assertTrue(a._testPipeFifoNoBlocking())
    # only if root in PATH, with $ROOTSYS
    if os.getenv("ROOTSYS") is None:
      logger.warning("no env var $ROOTSYS, no _testPipeFifoRootNoBlocking root test")
    else:
      self.assertTrue(a._testPipeFifoRootNoBlocking())

  def test_030(self):
    if verbose: logger.info("test_030")
    a = PFIFO.PipeFifo(_aPipeName)
    self.assertEqual(a.readNoBlocking(), "") # nothing receive

  def test_040(self):
    if verbose: logger.info("test_040")
    a = PFIFO.PipeFifo(_aPipeName)
    self.assertTrue(a.writeNoBlocking("Yeah 1", Verbose=verbose)) # avoid log pipe needs to be somewhere previously open read
    self.assertEqual(a.readNoBlocking(), "Yeah 1")
    self.assertTrue(a.writeNoBlocking("Yeah 2"))
    self.assertTrue(a.writeNoBlocking("Yeah 3"))
    self.assertEqual(a.readNoBlocking(), "Yeah 2Yeah 3")
    # multiple open read at same time accepted ?
    b = PFIFO.PipeFifo(_aPipeName)
    self.assertTrue(b.writeNoBlocking("Yooo 1"))
    self.assertEqual(b.readNoBlocking(), "Yooo 1")
    self.assertTrue(b.writeNoBlocking("Yooo 2"))
    self.assertTrue(b.writeNoBlocking("Yooo 3"))
    self.assertEqual(b.readNoBlocking(), "Yooo 2Yooo 3")

  def test_041(self):
    if verbose: logger.info("test_041")
    Verbose = False
    th2 = PFIFO.ThreadPipeFifoSend(None, _aPipeName)
    th2.verbose = verbose
    if Verbose: PFIFO.printListOfThreadsAlive()
    self.assertEqual(len(PFIFO.getListOfThreadsAlive()), 1)
    th2.start()
    if Verbose: PFIFO.printListOfThreadsAlive()
    time.sleep(2)
    self.assertEqual(len(PFIFO.getListOfThreadsAlive()), 2)
    time.sleep(2)
    th2.stop()
    th2.join()
    self.assertEqual(len(PFIFO.getListOfThreadsAlive()), 1)
    self.assertEqual(th2.nbSend, 0)
    self.assertEqual(th2.nbReceive, 0)
    if Verbose: PFIFO.printListOfThreadsAlive()

  def test_042(self):
    if verbose: logger.info("test_042")
    Verbose = False
    th1 = PFIFO.ThreadPipeFifoReceiveBlocking(None, _aPipeName)
    th2 = PFIFO.ThreadPipeFifoSend(None, _aPipeName)
    th1.verbose = verbose
    th2.verbose = verbose
    if Verbose: PFIFO.printListOfThreadsAlive()
    self.assertEqual(len(PFIFO.getListOfThreadsAlive()), 1)
    th2.start()
    th1.start()
    #th2.start()
    if Verbose: PFIFO.printListOfThreadsAlive()
    time.sleep(1)
    th2.sendMessage("ONE at " + th1.getStrTime())
    time.sleep(2)
    th2.sendMessage("TWO at " + th2.getStrTime())
    self.assertEqual(len(PFIFO.getListOfThreadsAlive()), 3)
    time.sleep(5)
    th1.stop()
    th2.stop()
    time.sleep(2) #thx have time to stop....
    self.assertEqual(th2.nbSend, 2)
    self.assertEqual(th1.nbReceive, 2)
    if Verbose: PFIFO.printListOfThreadsAlive()
    self.assertEqual(len(PFIFO.getListOfThreadsAlive()), 1)
    # no join() useful
    th1.join()
    th2.join()
    if Verbose: PFIFO.printListOfThreadsAlive()
    self.assertEqual(len(PFIFO.getListOfThreadsAlive()), 1)

  def test_044(self):
    if verbose: logger.info("test_044")
    Verbose = False
    th1 = PFIFO.ThreadPipeFifoReceiveBlocking(None, _aPipeName)
    th2 = PFIFO.ThreadPipeFifoSend(None, _aPipeName)
    th1.verbose = verbose
    th2.verbose = verbose
    if Verbose: PFIFO.printListOfThreadsAlive()
    self.assertEqual(len(PFIFO.getListOfThreadsAlive()), 1)
    #th2.start()
    th1.start()
    th2.start()
    if Verbose: PFIFO.printListOfThreadsAlive()
    time.sleep(5)
    th2.sendMessage("ONE at " + th1.getStrTime())
    time.sleep(2)
    th2.sendMessage("TWO at " + th2.getStrTime())
    self.assertEqual(len(PFIFO.getListOfThreadsAlive()), 3)
    time.sleep(2)
    th1.stop()
    th2.stop()
    time.sleep(2) #thx have time to stop....
    self.assertEqual(th2.nbSend, 2)
    self.assertEqual(th1.nbReceive, 2)
    if Verbose: PFIFO.printListOfThreadsAlive()
    self.assertEqual(len(PFIFO.getListOfThreadsAlive()), 1)
    # no join() useful
    th1.join()
    th2.join()
    if Verbose: PFIFO.printListOfThreadsAlive()
    self.assertEqual(len(PFIFO.getListOfThreadsAlive()), 1)

  def test_046(self):
    #test remove named pipe during threads
    if verbose:
      logger.info("test_046 TODO there is problem th2 continue...?, fix it")
    Verbose = False
    th1 = PFIFO.ThreadPipeFifoReceiveNoBlocking(None, _aPipeName)
    th1.verbose = False
    th2 = PFIFO.ThreadPipeFifoSend(None, _aPipeName)
    th1.className = "receive_046"
    th2.className = "send_046"
    th1.verbose = verbose
    th2.verbose = verbose
    if Verbose: PFIFO.printListOfThreadsAlive()
    self.assertEqual(len(PFIFO.getListOfThreadsAlive()), 1)
    #th2.start()
    th1.start()
    th2.start()
    if Verbose: PFIFO.printListOfThreadsAlive()
    time.sleep(2)
    th2.sendMessage("ONE at " + th1.getStrTime())
    time.sleep(2)
    self.assertEqual(len(PFIFO.getListOfThreadsAlive()), 3)
    self.removePipe() #remove named pipe file
    time.sleep(3) #thx have time to stop....
    self.assertEqual(th2.nbSend, 1)
    self.assertEqual(th1.nbReceive, 1)
    if Verbose: PFIFO.printListOfThreadsAlive()
    #th1.stop() # ? th2 continue...? stop it
    th2.stop() # ? th2 continue...? stop it
    time.sleep(3) #thx have time to stop....
    self.assertEqual(len(PFIFO.getListOfThreadsAlive()), 1)
    # no join() useful


if __name__ == '__main__':
  verbose = False #False #True #False #True
  #PFIFO.verbose = verbose # True #
  unittest.main()
  pass
