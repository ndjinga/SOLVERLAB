#!/usr/bin/env python
# -*- coding: utf-8 -*-

# %% LICENSE_SALOME_CEA_BEGIN
# see PACKAGESPY/LICENCE file
# %% LICENSE_END


import unittest
import sys
import xyzpy.loggingXyz as LOG

isMain = False

class TestCase(unittest.TestCase):
  def test_000(self):
    logger = LOG.getLogger()
    self.assertNotEqual(logger, None)
    #print dir(logger),logger.name
    self.assertTrue('Logger' in logger.__class__.__name__) # 'RootLogger' 'Logger'
    # self.assertTrue('StreamHandler' in LOG.consoleHandler.__class__.__name__) # 'StreamHandlerIra' etc.

  def obsolete_test_010(self):
    logger = LOG.getLogger()
    level = LOG.getLoggerLevelName()
    LOG.pushLevel(123)
    self.assertTrue('123' in LOG.getLoggerLevelName()) # 'Level 123'
    LOG.popLevel()
    self.assertEqual(LOG.getLoggerLevelName(), level)

    if isMain:
      #one pop without push do not raise error ... but risky avoid it if not isMain
      LOG.popLevel(warning=False)
      #depends of push/sets before (AllTestLauncher)...
      self.assertEqual(LOG.getLoggerLevelName(), "WARNING")

    LOG.pushLevel(LOG.levels["CRITICAL"]+1)
    self.assertTrue('51' in LOG.getLoggerLevelName()) # 'Level 51'
    LOG.popLevel()

  def obsolete_test_020(self):
    debug =False
    logger = LOG.getLogger()
    level = LOG.getLoggerLevelName()
    try:
      lev = int(level)
    except:
      self.assertEqual(level in list(LOG.levels.keys()), True)

    for i in list(LOG.levels.keys()):
      LOG.pushLevel(i)
      self.assertEqual(LOG.getLoggerLevelName() , i)
      LOG.popLevel()
      self.assertEqual(LOG.getLoggerLevelName() , level)

    if debug: print("original pushLevel %s" % (LOG.getLoggerLevelName()))
    LOG.pushLevel(123)
    if debug: print("first    pushLevel %s" % (LOG.getLoggerLevelName()))
    levels = list(LOG.levels.keys())
    for i in levels:
      LOG.pushLevel(i)
      if debug: print("pushLevel %s -> %s" % (i, LOG.getLoggerLevelName()))
      self.assertEqual(LOG.getLoggerLevelName() , i)
    for i in range(len(levels)):
      LOG.popLevel()
      if debug: print("popLevel %s" % LOG.getLoggerLevelName())
    self.assertTrue('123' in LOG.getLoggerLevelName())
    LOG.popLevel()
    self.assertEqual(LOG.getLoggerLevelName() , level)

if __name__ == '__main__':
  isMain = True
  unittest.main()
  pass
