#!/usr/bin/env python
# -*- coding: utf-8 -*-

# %% LICENSE_SALOME_CEA_BEGIN
# see PACKAGESPY/LICENCE file
# %% LICENSE_END


"""\
tests for incrementOnNamesXyz
"""

import os
import unittest
import xyzpy.incrementOnNamesXyz as IONA

verbose = False

class TestCase(unittest.TestCase):
  
  class TestXX(object):
    pass

  class TestYY(object):
    pass

  def test_010(self):
    a = self.TestXX()
    self.assertEqual(IONA.getIncrement(a), 'TestXX_0')
    self.assertEqual(IONA.getIncrement(a), 'TestXX_1')
    IONA.resetIncrement(a)
    self.assertEqual(IONA.getIncrement(a), 'TestXX_0')
    self.assertEqual(IONA.getIncrement(a), 'TestXX_1')
    
  def test_020(self):
    a = self.TestXX()
    self.assertEqual(IONA.getIncrement(a), 'TestXX_2')
    self.assertEqual(IONA.getIncrement(a), 'TestXX_3')
    IONA.resetIncrement(a)
    self.assertEqual(IONA.getIncrement(a), 'TestXX_0')
    self.assertEqual(IONA.getIncrement(a), 'TestXX_1')
    
  def test_030(self):
    a = self.TestYY()
    self.assertEqual(IONA.getIncrement(a), 'TestYY_0')
    IONA.resetIncrement(self.TestXX())
    self.assertEqual(IONA.getIncrement(a), 'TestYY_1')
    self.assertEqual(IONA.getIncrement(self.TestXX()), 'TestXX_0')
    
  def test_040(self):
    a = self.TestXX()
    b = self.TestYY()
    self.assertEqual(IONA.getIncrement(a), 'TestXX_1')
    self.assertEqual(IONA.getIncrement(b), 'TestYY_2')
    IONA.resetAllIncrements()
    self.assertEqual(IONA.getIncrement(a), 'TestXX_0')
    self.assertEqual(IONA.getIncrement(b), 'TestYY_0')
    
  def test_050(self):
    a = self.TestXX()
    b = self.TestYY()
    IONA.resetAllIncrements()
    self.assertEqual(IONA.getIncrement(a, "AA"), 'AA0')
    self.assertEqual(IONA.getIncrement(b, "BB"), 'BB0')
    
    self.assertEqual(IONA.getIncrement(a), 'AA1')
    self.assertEqual(IONA.getIncrement(b), 'BB1')
    
    self.assertEqual(IONA.getIncrement(a), 'AA2')
    self.assertEqual(IONA.getIncrement(b), 'BB2')
    
    self.assertEqual(IONA.getIncrement(a, "XXAA"), 'XXAA3')
    self.assertEqual(IONA.getIncrement(b, "XXBB"), 'XXBB3')

  def test_050(self):
    a = self.TestXX()
    b = self.TestYY()
    IONA.resetAllIncrements()
    #AffectedNames are yet from geom names from hdf5 (for example...)
    IONA.setAffectedNames(a, ['AA0', 'AA2'])
    IONA.setAffectedNames(b, ['BB1', 'BB3'])
    self.assertEqual(IONA.getIncrement(a, "AA"), 'AA1')
    self.assertEqual(IONA.getIncrement(b, "BB"), 'BB0')
    self.assertEqual(IONA.getIncrement(a, "AA"), 'AA3')
    self.assertEqual(IONA.getIncrement(b, "BB"), 'BB2')
    self.assertEqual(IONA.getIncrement(a, "AA"), 'AA4')
    self.assertEqual(IONA.getIncrement(b, "BB"), 'BB4')

if __name__ == '__main__':
  unittest.main()


