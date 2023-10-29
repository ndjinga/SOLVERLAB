#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import unittest
import objargpy.objarg as OARG

verbose = False
   
#########################################
class TestCase(unittest.TestCase):

  def test_005(self):
    a = OARG.ObjArg()
    if verbose: a.printHelp()
    self.assertTrue("show this help message and exit" in a.getHelp())
    self.assertTrue("--perfectSquare" in a.getHelp())
 
  def test_010(self):
    a = OARG.ObjArg()
    if verbose: print("test_010\n%s" % a)
    self.assertTrue('2.2' in str(a))
    self.assertTrue('[123, 456.7, 8910]' in str(a))

  def test_015(self):
    a = OARG.ObjArg("-v -ps 16".split())
    if verbose: print("test_015\n%s" % a)
    self.assertEqual(a.perfectSquare, 16)
    self.assertEqual(a.perfectSquare, 16)
    a.perfectSquare = 11
    a.getInitialNamespace().perfectSquare = 11
    self.assertEqual(a.perfectSquare, 11)
    self.assertEqual(a.getInitialNamespace().perfectSquare, 16)
    with self.assertRaises(Exception):
      OARG.ObjArg("-v -ps 11".split())

  def test_020(self):
    with self.assertRaises(Exception):
      OARG.ObjArg("-v -f 1.2".split())
    with self.assertRaises(Exception):
      OARG.ObjArg("-v -f 34".split())
    with self.assertRaises(Exception):
      OARG.ObjArg("-v -ps 11".split())
    with self.assertRaises(Exception):
      OARG.ObjArg("-unknown 34".split())

  def test_025(self):
    a = OARG.ObjArg()
    b = OARG.ObjArg()
    #"-t" or "--thirdAttList"
    a.bbb = OARG.ObjArg("--t 11,22,44".split())
    if verbose: print("test_025\n%s" % str(a))
    self.assertTrue('[123, 456.7, 8910]' in str(a))
    self.assertTrue('[123, 456.7, 8910]' in str(a.thirdAttList))
    self.assertFalse('[11, 22, 44]' in str(a))
    if verbose: print("test_025\n%s" % a.reprAll())
    self.assertTrue('[123, 456.7, 8910]' in a.reprAll())
    self.assertTrue('[11, 22, 44]' in a.reprAll())
    self.assertTrue('[11, 22, 44]' in str(a.bbb))
    #a.reprAll(), is something like res, here for info
    res = """
ObjArg([
  ('bbb', ObjArg([
    ('firstAttInt', 1),
    ('perfectSquare', 9),
    ('secondAttFloat', 2.2),
    ('thirdAttList', [11, 22, 44]),
    ('verbose', False),
  ]),
  ('firstAttInt', 1),
  ('perfectSquare', 9),
  ('secondAttFloat', 2.2),
  ('thirdAttList', [123, 456.7, 8910]),
  ('verbose', False),
])"""

    
#########################################
# main
#########################################

if __name__ == '__main__':
  if verbose: print("objargTest.py: sys.argv:", sys.argv)
  verbose = True
  unittest.main()
  pass
   

