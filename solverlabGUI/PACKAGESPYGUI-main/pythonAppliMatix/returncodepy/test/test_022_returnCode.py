#!/usr/bin/env python
#-*- coding:utf-8 -*-

# %% LICENSE_SALOME_CEA_BEGIN
# see PACKAGESPY/LICENCE file
# %% LICENSE_END



import os
import sys
import unittest
import pprint as PP

import debogpy.debug as DBG
from returncodepy.returnCode import ReturnCode as RC

verbose = False # True

class TestCase(unittest.TestCase):
  "Test the debug.py"""
  
  def test_000(self):
    # one shot setUp() for this TestCase
    if verbose:
      DBG.push_debug(True)
      DBG.write("assert unittest", [a for a in dir(self) if "assert" in a])
    pass
  
  def test_010(self):
    rc = RC()
    self.assertFalse(rc.isOk())
    rrc = str(rc)
    DBG.write("test_010 str", rrc)
    self.assertIn("ND:", rrc)
    self.assertIn("No given explanation", rrc)
    self.assertNotIn("--value", rrc)
    rrc = repr(rc)
    DBG.write("test_010 repr", rrc)
    self.assertIn("ND:", rrc)
    self.assertIn("No given explanation", rrc)
    self.assertIn("--value", rrc)
       
  def test_015(self):
    rc = RC("OK", "all is good")
    self.assertTrue(rc.isOk())
    rrc = str(rc)
    DBG.write("test_015 str", rrc)
    self.assertIn("OK:", rrc)
    self.assertIn("all is good", rrc)
    self.assertNotIn("--value", rrc)
    rrc = repr(rc)
    DBG.write("test_015 repr", rrc)
    self.assertIn("OK:", rrc)
    self.assertIn("all is good", rrc)
    self.assertIn("None", rrc)
    aVal = "I am a value result"
    rc.setValue(aVal)
    self.assertTrue(rc.isOk())
    self.assertEqual(rc.getValue(), aVal)
    rrc = repr(rc)
    DBG.write("repr", rrc)
    
  def test_020(self):
    aVal = "I am a value result"
    rc1 = RC("OK", "all is good1", aVal + "1")
    self.assertTrue(rc1.isOk())
    rc2 = RC("OK", "all is good2", aVal + "2")
    self.assertTrue(rc2.isOk())
    rc3 = rc1 + rc2
    self.assertTrue(rc3.isOk())
    rrc = repr(rc3)
    DBG.write("test_020 repr", rrc)
    self.assertIn("OK:", rrc)
    self.assertIn("good1", rrc)
    self.assertIn("good2", rrc)
    self.assertIn("result1", rrc)
    self.assertIn("result2", rrc)
    rc4 = rc3 + rc1
    rrc = repr(rc4)
    DBG.write("test_020 repr", rrc)
    self.assertEqual(len(rc4.getWhy()), 3)
    self.assertEqual(len(rc4.getValue()), 3)
    
  def test_025(self):
    rc0 = RC("KO")
    aVal = "I am a value result"
    rc1 = RC("OK", "all is good1", aVal + "1")
    self.assertTrue(rc1.isOk())
    rc1.setStatus("KO") # raz status and why and value
    self.assertFalse(rc1.isOk())
    self.assertEqual(repr(rc0), repr(rc1))
    
  def test_026(self):
    rc0 = RC("KO")
    aVal = "I am a value result"
    rc1 = RC("OK", "all is good1", aVal)
    rc2 = rc0 + rc1 # make list with two why and two value
    DBG.write("test_026 str", str(rc2))
    DBG.write("test_026 repr", repr(rc2))
    self.assertFalse(rc2.isOk())
    self.assertIn("KO:", repr(rc2))
    self.assertIn("KO:", str(rc2))    
    self.assertEqual(len(rc2.getWhy()), 2)
    self.assertEqual(len(rc2.getValue()), 2)

  def test_027(self):
    rc0 = RC("KO")
    aVal = "I am a value result"
    rc1 = RC("OK", "all is good1", aVal)
    rc2 = sum([rc0, rc1]) # make list with two why and two value
    self.assertFalse(rc2.isOk())
    self.assertIn("KO:", repr(rc2))
    self.assertIn("KO:", str(rc2))    
    self.assertEqual(len(rc2.getWhy()), 2)
    self.assertEqual(len(rc2.getValue()), 2)
    
  def test_028(self):
    rc0 = RC("OK")
    aVal = ["I am a mutable list value result"]
    rc1 = RC("OK", "all is good1", aVal)
    rc2 = rc0 + rc1 + rc1 + rc0 # make list with why and value
    DBG.write("test_028 repr", str(rc2))
    DBG.write("test_028 repr", repr(rc2))
    self.assertTrue(rc2.isOk())
    self.assertIn("OK:", repr(rc2))
    self.assertIn("OK:", str(rc2))    
    self.assertEqual(len(rc2.getWhy()), 4)
    self.assertEqual(len(rc2.getValue()), 4)
    
    aVal[0] = "modified mutable list value result"
    DBG.write("test_028 repr", repr(rc1))
    DBG.write("test_028 repr", repr(rc2))
    # deepcopy value no clearly assumed, could be tricky
    self.assertIn("modified mutable", repr(rc1))
    self.assertNotIn("modified mutable", repr(rc2))

          
  def test_999(self):
    # one shot tearDown() for this TestCase
    if verbose:
      DBG.pop_debug()
    return
    
if __name__ == '__main__':
    unittest.main(exit=False)
    pass

