#!/usr/bin/env python
#-*- coding:utf-8 -*-

# %% LICENSE_SALOME_CEA_BEGIN
# see PACKAGESPY/LICENCE file
# %% LICENSE_END


import os
import sys
import platform
import unittest

import test.initializeTest # set PATH etc for test


class TestCase(unittest.TestCase):
  "Test the debug.py"""
  
  def test_000(self):
    # one shot setUp() for this TestCase
    # DBG.push_debug(True)
    # SAT.setNotLocale() # test english
    return
    
  def test_005(self):
    import debogpy.debug as DBG # Easy print stderr (for DEBUG only)
    res = DBG.getLocalEnv()
    self.assertTrue(len(res.split()) > 0)
    self.assertTrue("USERNAME :" in res)
    if platform.system() != "Windows":
      self.assertTrue("LANG :" in res)
       
  # TODO another test without pyconf
      
  def test_999(self):
    # one shot tearDown() for this TestCase
    # SAT.setLocale() # end test english
    # DBG.pop_debug()
    return
    
if __name__ == '__main__':
    unittest.main(exit=False)
    pass

