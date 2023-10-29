#!/usr/bin/env python
#-*- coding:utf-8 -*-

#  Copyright (C) 2010-2023  CEA
#
#  This library is free software; you can redistribute it and/or
#  modify it under the terms of the GNU Lesser General Public
#  License as published by the Free Software Foundation; either
#  version 2.1 of the License.
#
#  This library is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#  Lesser General Public License for more details.
#
#  You should have received a copy of the GNU Lesser General Public
#  License along with this library; if not, write to the Free Software
#  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA


import os
import sys
import unittest
import pprint as PP

import debogpy.debug as DBG
import solverlabpy.dateTime as DATT
import time as TI

verbose = False #True

class TestCase(unittest.TestCase):
  "Test the debug.py"""
  
  def test_000(self):
    # one shot setUp() for this TestCase
    if verbose:
      DBG.push_debug(True)
      # DBG.write("assert unittest", [a for a in dir(self) if "assert" in a])
    pass
  
  def test_010(self):
    t = DATT.DateTime()
    self.assertFalse(t.isOk())
    rrt = str(t)
    DBG.write("test_010 str", rrt)
    self.assertIn("Undefined", rrt)
    rrt = repr(t)
    DBG.write("test_010 repr", rrt)
    self.assertIn("DateTime", rrt)
    self.assertIn("Undefined", rrt)
    with self.assertRaises(Exception):
      t.raiseIfKo()

       
  def test_015(self):
    t = DATT.DateTime("now")
    self.assertTrue(t.isOk())
    rrt = str(t)
    DBG.write("test_015 str", rrt)
    self.assertIn("20", rrt) # 2018 to 2099 ok
    self.assertIn("-", rrt)
    self.assertIn(":", rrt)
    rrt = repr(t)
    DBG.write("test_015 repr", rrt)
    self.assertIn("DateTime", rrt)
    self.assertIn("20", rrt) # 2018 to 2099 ok
    self.assertIn("-", rrt)
    self.assertIn(":", rrt)

    
  def test_020(self):
    t1 = DATT.DateTime("now")
    t2 = DATT.DateTime(t1)
    self.assertTrue(t2.isOk())
    self.assertEqual(t1, t2)
    t2 = DATT.DateTime("now")
    self.assertNotEqual(t1, t2) # microseconds differs
    
    DATT.sleep(3) # 3 second more
    t2 = DATT.DateTime("now")
    self.assertGreater(2, 1) # to be sure
    self.assertGreater(str(t2), str(t1)) # seconds differs
    self.assertGreater(repr(t2), repr(t1)) # seconds differs
    self.assertGreater(t2, t1)
    self.assertTrue(t2 > t1)
    self.assertFalse(t2 == t1)
    self.assertFalse(t2 < t1)
    self.assertFalse(t2 <= t1)
    
  def test_040(self):
    t1 = DATT.DateTime("now")
    delta = DATT.DeltaTime(t1)
    self.assertFalse(delta.isOk())
    self.assertIn("Undefined", delta.toSeconds()) 
    DBG.write("test_040 str", str(delta))
    DBG.write("test_040 repr", repr(delta))   
    with self.assertRaises(Exception):
      delta.raiseIfKo()
      DATT.DateTime().raiseIfKo()
       
  def test_042(self):
    t1 = DATT.DateTime("now")
    DATT.sleep(3.1) # 3.1 second more
    t2 = DATT.DateTime("now")
    self.assertTrue(t2 > t1)
    delta = DATT.DeltaTime(t1, t2)
    self.assertGreater(delta.toSeconds(), 3)
    self.assertEqual(int(delta.toSeconds()), 3)
    DBG.write("test_042 str", str(delta))
    DBG.write("test_042 repr", repr(delta))
    delta2 = delta.raiseIfKo()
    self.assertEqual(delta2.toSeconds(), delta.toSeconds())

  def test_044(self):
    t1 = DATT.DateTime("now")
    t2 = DATT.DateTime(t1) + 3.1 # add 3 seconds
    delta = DATT.DeltaTime(t1, t2)
    self.assertGreater(delta.toSeconds(), 3)
  
    ti = TI.time()
    t4 = DATT.DateTime(ti)
    t5 = t4 + 10.1
    DBG.write("test_044 ti", [type(ti), ti, t4, t5])
    delta = DATT.DeltaTime(t4, t5)
    self.assertGreater(delta.toSeconds(), 10)
    DBG.write("test_044 delta", [delta.toStrHuman(), delta.toStrHms()])
    self.assertNotIn("-", delta.toStrHuman())
    self.assertNotIn("-", delta.toStrHms())
    
    delta = DATT.DeltaTime(t5, t4) # negative delta
    self.assertLess(delta.toSeconds(), -10)
    DBG.write("test_044 delta", [delta.toStrHuman(), delta.toStrHms()])
    self.assertIn("-10s", delta.toStrHuman())
    self.assertIn("-0h0m10s", delta.toStrHms())
    
  def test_046(self):
    for more in [0, 0.56789, 5.6789, 56.789, 61, 3661, 36061]:
      t1 = DATT.DateTime("now")
      t2 = DATT.DateTime(t1)
      t2.addSeconds(more)
      delta = DATT.DeltaTime(t1, t2)
      r = delta.toStrHuman()
      DBG.write("test_046 str", r)
      if more < 60: 
        self.assertIn("s", r)
        self.assertNotIn("m", r)
        self.assertNotIn("h", r)
        continue
      if more < 3600: 
        self.assertIn("s", r)
        self.assertIn("m", r)
        self.assertNotIn("h", r)
      else:
        self.assertIn("s", r)
        self.assertIn("m", r)
        self.assertIn("h", r)

      
          
  def test_999(self):
    # one shot tearDown() for this TestCase
    if verbose:
      DBG.pop_debug()
    return
    
if __name__ == '__main__':
    unittest.main(exit=False)
    pass

