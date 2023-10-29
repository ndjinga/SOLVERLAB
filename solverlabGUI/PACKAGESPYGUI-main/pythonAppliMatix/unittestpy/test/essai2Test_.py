#!/usr/bin/env python
#-*- coding:utf-8 -*-

"""
This test creates error explicitely for test html unitest etc,
so it is named  xxx_.test to be run explicitely
"""

import unittest
import HTMLTestRunner as HTST

verbose = False

class TestCase(unittest.TestCase):

  def test_005(self):
    """essai2 test_005 for doc"""
    if verbose: print('hello test_005')
    a="  123456  "
    self.assertEqual(a.strip(), "Ooops")
   
  def test_010(self):
    """essai2 test_010 for doc"""
    if verbose: print('hello test_010')
    a="  123   456  "
    self.assertNotEqual(a, "123456")
   
  def test_015(self):
    """essai2 test_015 for doc"""
    if verbose: print('hello test_010')
    OOps22222
   
   
if __name__ == '__main__':
  verbose = True
  #unittest.main(verbosity=True)
  HTST.main()
  pass
