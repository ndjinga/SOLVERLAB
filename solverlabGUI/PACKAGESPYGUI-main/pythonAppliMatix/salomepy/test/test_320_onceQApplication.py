#!/usr/bin/env python
# -*- coding: utf-8 -*-

# %% LICENSE_SALOME_CEA_BEGIN
# see PACKAGESPY/LICENCE file
# %% LICENSE_END


import unittest
import sys
from salomepy.onceQApplication import OnceQApplication
import xyzpy.loggingXyz as LOG

class TestCase(unittest.TestCase):
  def test_000(self):
    LOG.pushLevel("CRITICAL")
    app1=OnceQApplication(sys.argv)
    self.assertNotEqual(app1, None)
    self.assertEqual(app1.__class__.__name__, 'QApplication')
    app2=OnceQApplication(sys.argv)
    self.assertEqual(app1, app2)
    self.assertEqual(id(app1), id(app2))
    app3=OnceQApplication()
    self.assertEqual(app1, app3)
    self.assertEqual(id(app1), id(app3))
    LOG.popLevel()

if __name__ == '__main__':
  """OnceQApplicationTest.py"""
  unittest.main()
  pass

